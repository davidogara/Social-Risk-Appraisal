; ********************************************************************************************************************
; Model: Social Risk Appraisal
; team: David O'Gara, Matt Kasman, Laurent Hebert-Dufrense, and Ross A. Hammond
; author: David O'Gara | Division of Computational and Data Sciences | Washington University in St. Louis
; last updated: 2024-01-30
; ********************************************************************************************************************


;extensions [vid]
extensions [py nw table]
globals [count-disease-s count-disease-i count-disease-r percent-infected percent-hiding percent-fearful
  _recording-save-file-name color-min color-max count-fearful count-hiding day-of-max-infected count-disease-I-prev gen-random-number random-count
  infinity apl-nw cc-nw
sort-idx-to-turtle-idx
turtle-idx-to-sort-idx
original-sort-idx-to-turtle-idx
fear-score-of-infected
fear-score-of-infected-cumulative
avg-hide-time
]

turtles-own [

  disease-state ;; S, I, or R
  fear-state ;; S, I, or R
  time-infected
  time-recovering
  time-fearful
  hiding?
  time-hidden
  fear-count ; for complex contagion
  fear-score ; for FTA
  fear-ever-infected?
  disease-rate ; Disease-Infection-Rate
  tick-fearful ; tick at which an agent became infected
  anchor? ; are they an influential agent?
  fear-score-step-size
  my-clustering-coefficient ; for clustering
  distance-from-other-turtles
  fear-score-sort
  time-to-hide ; only used for FTA-v2
  hide-score
  initial-fear-score ; record fear-score upon model instantiation
  pr-release
  can-spread-fear?
  fear-score-when-infected
  tick-infected
  dend ; dendogram to trace infection
]

links-own [
  rewired?
]

to init-colors
  set color-min min [fear-score] of turtles
  set color-max max [fear-score] of turtles
  ask patches [set pcolor 99]
end

to color-patches-based-on-fear-score
  let denom (color-max - color-min)
  ;set denom max (list denom 1)
  ask turtles [
    let my-fear fear-score
    let newcolor 99 - (9 * (my-fear / denom))
    ask neighbors [
      set pcolor newcolor
    ]
    ask patch-here [set pcolor newcolor]
  ]

end


to setup

  clear-all
  reset-ticks
  if( Random-Seed? = true)[
    random-seed behaviorspace-run-number + random_seed_num
  ]
  ;;Set Globals
  set count-disease-I-prev 0
  set infinity 9999 ; arbitrary large number

  ask n-of (num-turtles) patches
  [ sprout 1
    [ set disease-state "S"
      set fear-state "S"
      set hiding? false
      set color blue
      set shape "person"
      set time-infected 0
      set time-fearful 0
      set time-recovering 0
      set time-hidden 0
      set time-to-hide 0
      set anchor? false
      set fear-score-step-size step-increment
      set can-spread-fear? true
      set tick-infected -100 ; to account for agents being infected at tick 0
      set dend "No-infections-yet"
    ]
  ]

  ; why do we do this? Because if we ask init-infected patches to sprout 1, you cluster the initial infecteds when forming the `Small World`
  ask n-of (init-infected) turtles
  ;ask one-of patches
    [ set disease-state "I"
      set fear-state "S"
      set hiding? false
      set color red
      set shape "person"
      set time-infected 0
      set time-fearful 0
      set time-recovering 0
      set time-hidden 0
      set time-to-hide 0
      set anchor? false
      set fear-score-step-size step-increment
    ]
  ;; Fear-Prob setup ;;
  if (type-of-fear-distribution = "Gamma") [setup-fear-mechanism]
  if (type-of-fear-distribution = "Anchor")[create-trunc-normal-with-anchors]
  if (type-of-fear-distribution = "Truncated Normal") [create-trunc-normal]
  if (type-of-fear-distribution = "Custom Normal") [create-custom-normal-with-middle-carved-out]
  if (type-of-fear-distribution = "Custom Normal v2") [create-custom-normal-separated-clusters]
  if (type-of-fear-distribution = "Uniform") [ask turtles [set fear-score random-float 1.0]]

  ask turtles [
    set initial-fear-score fear-score
    set fear-score-when-infected fear-score
  ]

  ;;Network Setup;;;
  if network-type = "Spatially Clustered" [setup-spatially-clustered-network]
  if network-type = "Random" [setup-random-network]
  if network-type = "Small World" [setup-small-world-network]

  if network-type = "Small World Fear-Score" [setup-homophily-small-world]

  ask turtles with [count link-neighbors = 0][
    let choice (min-one-of (other turtles with [not link-neighbor? myself])
                   [distance myself])
      if choice != nobody [ create-link-with choice ]
  ]


  init-colors
  setup-disease
  update-globals
  do-plots


  set avg-hide-time []
  if (export-network? = true and behaviorspace-experiment-name != "")[
  let filename (word "networks-for-assortativity/" behaviorspace-experiment-name "_" behaviorspace-run-number ".csv")
  export-world filename
  ]
end




to go
  ;move-turtles
  update-disease-states
  update-fear-states
  update-removal-from-hiding
  update-can-spread-fear
  ask turtles [update-colors]
  update-globals
  do-plots
  tick
  if (ticks >= 150) [
    stop
  ]
  if (color-patches?  = true)[color-patches-based-on-fear-score]
  reset-hiding
end


to setup-disease

  ask turtles [
    set disease-rate Disease-Infection-Rate
  ]

end

to setup-fear-mechanism
  if (Fear-Mechanism = "FTA" or Fear-Mechanism = "Complex Contagion-Gamma" or Fear-Mechanism = "FTA-Mean" or Fear-Mechanism = "FTA-v2") [
    ask turtles [
      set fear-score random-beta Fear-Infection-Rate (1 - Fear-Infection-Rate)
    ]
  ]

  if (Fear-Mechanism = "Complex Contagion")[
    if (Random-Seed? = true)[
      ; generate *two* random numbers to stay in sync with the random-seed under different fear mechanisms
      ask turtles [set fear-score random-beta Fear-Infection-Rate (1 - Fear-Infection-Rate)]
    ]
    ; and then set it to the thing we *actually* want
      ask turtles [set fear-score Fear-Infection-Rate]
  ]

  ask turtles [
    set fear-score fear-score / beta-dist-denom
  ]

end



to-report random-beta [ #alpha #beta ]
  let XX random-gamma #alpha 1
  let YY random-gamma #beta 1
  report XX / (XX + YY)
end

;; Disease-States
to update-disease-states

  ;;Spread Disease and Fear from Non-Hiding Disease-Infecteds. -- Being sick and hiding does not spread fear or disease
  ask-concurrent turtles with [disease-state = "I" and hiding? = false] [
    let infector who
    ask link-neighbors [
      ;;Disease Transmission
      if (disease-state = "S" and hiding? = false)[
        let prob random-float 1
        if (prob < disease-rate) [ ;;Disease-Infection-Rate)[
          set disease-state "I"
          set fear-score-when-infected fear-score
          set tick-infected ticks
          set dend (word "I-got-infected-by-" infector "-and-the-prob-was-" prob)
        ]
      ]
      ; Fear Transmission Through Diseased Neighbors
      ; Fear-Infecteds can spread fear via simple contagion, rather than complex, because seeing the "sick" is not the same
      ; as just seeing the "scared"
      if (Fear-On?)[
        if (fear-state = "S" and hiding? = false)[
          if (random-float 1 < fear-score)[
            ;let x 1
            set fear-state "I"
            set time-fearful 0
            set tick-fearful ticks
          ]
        ]
      ]
    ]
  ]

  ;;Recover/Update All Disease-Infecteds.
  ask-concurrent turtles with [disease-state = "I"] [
    ifelse (time-infected >= infection-duration)
    [
      set disease-state "R"
    ]
    [
      set time-infected time-infected + 1
    ]
  ]
end


;; Fear States;;

to update-fear-states
  update-fear-states-general
end




to update-fear-states-general
  if (Fear-On?)[
    ;;Spread Fear from Fear-Infecteds (only those who became afraid before this round).
    ask-concurrent turtles with [fear-state = "I" and time-fearful > 0 and can-spread-fear?][
      ask link-neighbors [
        ;;Fear Transmission Through Fearful Neighbors
        if (fear-state = "S" and count link-neighbors with [fear-state = "I"] >= fear-threshold)[
            if (random-float 1.0 < fear-score)[
               set fear-state "I"
               set time-fearful 0
               set time-hidden 0
               set hiding? True
               set fear-ever-infected? True
               set tick-fearful ticks
          ]

          ]

      ]
    ]

    ;;Recover/Update All Fear-Infecteds
    ask turtles with [fear-state = "I"] [

      ifelse (time-fearful >= hide-duration)
      [
        set fear-state "S"
        set time-fearful 0
      ]
      [
        set time-fearful time-fearful + 1
      ]
    ]
  ]



  ;; Update Fear-Score if Mechanism is FTA

  if (Fear-Mechanism = "FTA") [ask turtles [change]]
end


;; Fear-Score Update rule
to change

  let radius 0
  let increment 0
  set increment min sentence fear-score-step-size (abs (mean [fear-score] of link-neighbors - fear-score))

  if(satisficing?)[
    set radius satisficing-radius
  ]

  ;; Allowing for errors,
  carefully [

    ;; Check whether neighbors' mean fear-score differs perceptibly from yours. If so, adjust
    ;;   fear-score
    ifelse (fear-score < mean[fear-score - radius] of link-neighbors) [set fear-score fear-score + increment]
    [
      if fear-score > mean[fear-score + radius] of link-neighbors [set fear-score fear-score - increment]
    ]

  ][let x 1];print error-message]

end

; helper function to make tanh
to-report tanh [x]
  let num exp (2 * x) - 1
  let denom exp (2 * x) + 1
  report num / denom

end

to update-removal-from-hiding
  ;;Increment Time-Hidden / Come Out of Hiding
  if (Release-Mechanism = "fixed-duration")[
  ask turtles with [hiding?] [
    ifelse (time-hidden >= hide-duration)
    [
      set hiding? false
      set fear-state "S"
      set time-hidden 0
      set time-fearful 0
    ]
    [
      set time-hidden time-hidden + 1
    ]
  ]
  ]
  if (Release-Mechanism = "probabilistic")[
  ;;Increment Time-Hidden / Come Out of Hiding
  ask turtles with [hiding?] [

    ; use a function of the form tanh
    let x (1 / (prob-release-damping * fear-score)) * time-hidden
    ;let prob-release tanh x

    ; use exponential decay
    let prob-hide exp(-1.0 * x)
    ;set pr-release prob-release

    ifelse(random-float 1.0 < prob-hide) [
       set time-hidden time-hidden + 1
       set time-fearful time-fearful + 1
    ][
      set hiding? false
      set fear-state "S"
      set avg-hide-time insert-item 0 avg-hide-time time-hidden
      set time-hidden 0
      set time-fearful 0
    ]
  ]
  ]

end


to reset-hiding
  ; let agents hide a again, if enough time has passed

  ask turtles with [(ticks - tick-fearful) >= restore-hiding and disease-state = "S" and fear-state = "S" and hiding? = false][
  set time-fearful 0
  set time-hidden 0
  ]

end

to update-can-spread-fear
  if (agents-can-spread-fear-while-hiding? = false)[
    ask turtles [
      if-else (hiding?)[
        set can-spread-fear? false
      ][
        set can-spread-fear? true
      ]
    ]
  ]

end

to update-colors
  if (disease-state = "S")
  [
    if (fear-state = "S")
    [set color blue]
    if (fear-state = "I")
    [set color yellow]
    if (fear-state = "R")
    [set color blue]
  ]
  if (disease-state = "I")
  [
    if (fear-state = "S")
    [set color red]
    if (fear-state = "I")
    [set color orange]
    if (fear-state = "R")
    [set color red]
  ]
  if (disease-state = "R")
  [
    set color white
  ]
  if (hiding?)
  [
    set color green
  ]

end


to update-globals
  ;; Counts
  set count-disease-S count (turtles with [disease-state = "S"])
  set count-disease-I count (turtles with [disease-state = "I"])
  set count-disease-R count (turtles with [disease-state = "R"])
  set count-fearful count (turtles with [fear-state = "I"])
  set count-hiding count (turtles with [hiding? = true])
  ;; Percents
  set percent-hiding count (turtles with [hiding? = true]) / num-turtles * 100
  set percent-fearful count (turtles with [fear-state = "I"]) / num-turtles * 100
  set percent-infected ((count-disease-I + count-disease-R)/ num-turtles * 100)

  ; Fear-Score of Infected
  set fear-score-of-infected 0
  if (count-disease-I > 0)[
    set fear-score-of-infected mean [initial-fear-score] of turtles with [disease-state = "I"]
  ]
  set fear-score-of-infected-cumulative mean [initial-fear-score] of turtles with [disease-state != "S"]

  ; update day of max-infected
  if (count-disease-I > count-disease-I-prev)[
    set day-of-max-infected ticks
    set count-disease-I-prev count-disease-I
  ]
end

to add-edge
  let node2 one-of turtles
  ifelse link-neighbor? node2 or node2 = self
    ;; if there's already an edge there, then go back
    ;; and pick new turtles
    [ add-edge ]
    ;; else, go ahead and make it
    [create-link-with node2]

end


to set-custom-fear-dist
  ask turtles [
        set fear-score random-normal 0.3 0.05
    if (random-float 1.0 < 0.5)[
        set fear-score random-normal 0.1 0.05
    ]
  ]

end



to do-plots
  if (behaviorspace-run-number = 0)[
  ; Plot #1
  set-current-plot "Disease SIR Curves"
  set-current-plot-pen "S"
  plot count-disease-S
  set-current-plot-pen "I"
  plot count-disease-I
  set-current-plot-pen "R"
  plot count-disease-R


  ; Plot #2
  set-current-plot "Number Hiding"
  plot count turtles with [hiding? = True]

  ; Plot #3
  set-current-plot "Count Fearful"
  plot count turtles with [fear-state = "I"]

  set-current-plot "Fear-Score Distribution"
  histogram [fear-score] of turtles

  ]

  ; Plot #4
  set-current-plot "Avg Fear Score of Infected"
  plot fear-score-of-infected-cumulative


end


to create-custom-normal-separated-clusters
  create-custom-normal-with-middle-carved-out
  let barrier-size 5
  ; move turtles
  ask turtles [
    if (fear-score < Fear-Infection-Rate)[
      move-to one-of patches with [(pxcor - pycor) < -1.0 * barrier-size]
    ]
    if (fear-score > Fear-Infection-Rate)[
      move-to one-of patches with [(pxcor - pycor) > barrier-size]
    ]


  ]

end


to create-custom-normal-with-middle-carved-out
  ;; make a custom normal distribution
  let window-size 0.05
  ask turtles [
  set fear-score random-normal Fear-Infection-Rate 0.05
  ; bound
  set fear-score min (list fear-score (2 * Fear-Infection-Rate))
  set fear-score max (list fear-score 0.0)
  ; remove middle
    if (fear-score < Fear-Infection-Rate and fear-score > (Fear-Infection-Rate - window-size))[
      set fear-score Fear-Infection-Rate - window-size
    ]
    if (fear-score > Fear-Infection-Rate and fear-score < (Fear-Infection-Rate + window-size))[
      set fear-score Fear-Infection-Rate + window-size
    ]

  ]


;  while [i < 50000][
;    let x random-normal 0.2 0.2 ; a normal dist with a little extra mass on each end
;    set x min (list x 0.4)
;    set x max (list x 0.0)
;    if
;
;    if ((x < 0.15 or x > 0.25) and (x > 0 and x < 0.4))[ ; for symmetry
;    set distribution fput x distribution
;      set i i + 1
;    ]
;  ]
;
;  ask turtles [
 ;   set fear-score one-of distribution
    ; if fear-score in tail [set fear-score-step-size 1e-6]
 ; ]
  let lower calc-pct 2.5 [fear-score] of turtles ;with [fear-score < 0.15]
  let upper calc-pct 97.5 [fear-score] of turtles ;with [fear-score > 0.25]

  ask turtles [
    if (fear-score <= lower or fear-score >= upper)[
      let tmp 1
      set anchor? true
      set fear-score-step-size 1e-6
    ]
  ]

  show lower
  show upper
end

to-report calc-pct [ #pct #vals ]
  let #listvals sort #vals
  let #pct-position #pct / 100 * length #vals
  ; find the ranks and values on either side of the desired percentile
  let #low-rank floor #pct-position
  let #low-val item #low-rank #listvals
  let #high-rank ceiling #pct-position
  let #high-val item #high-rank #listvals
  ; interpolate
  ifelse #high-rank = #low-rank
  [ report #low-val ]
  [ report #low-val + ((#pct-position - #low-rank) / (#high-rank - #low-rank)) * (#high-val - #low-val) ]
end

to create-separating-dist
  ; 0.1 and 0.3

  ask turtles [
    set fear-score random-normal Fear-Infection-Rate 0.05
    if (random-float 1.0 < 0.05)[
        set fear-score 0.001
        set anchor? true
        set fear-score-step-size 1e-6
      ]
    if (fear-score < 0)[
      set fear-score 0
    ]

  ]

end


to create-trunc-normal-with-anchors

  create-trunc-normal

  ask turtles with [fear-score = 0.01][
    if (random-float 1.0 < 0.25)[
      set anchor? true
      set fear-score-step-size 1e-6
    ]
  ]
end

to create-trunc-normal
  ; 0.1 and 0.3

  ask turtles [
    set fear-score random-normal Fear-Infection-Rate 0.05
    if (fear-score < 0.01)[
      set fear-score 0.01
    ]

  ]
  ask turtles [
  set fear-score fear-score / beta-dist-denom
  ]

  ; put a little mass out at ~0.5
  ask turtles [
    if (random-float 1.0 < 0.01)[
        ;set fear-score 0.5
        let x 1
      ]
  ]


end

;; NETWORK FUNCTIONS ;;

to setup-spatially-clustered-network
  let num-links (average-node-degree * num-turtles) / 2
  while [count links < num-links ]
  [
    ask one-of turtles
    [
      let choice nobody
      set choice (min-one-of (other turtles with [not link-neighbor? myself])
                   [distance myself])


      if (choice != nobody) [ create-link-with choice ]
    ]
  ]

  repeat 10
  [
    layout-spring turtles links 0.3 (world-width / (sqrt num-turtles)) 1
  ]
end

to setup-random-network
  let num-links (average-node-degree * num-turtles) / 2
  while [count links < num-links]
  [
    ask one-of turtles
    [
      let choice one-of turtles
      if choice != self [create-link-with choice]
    ]
  ]
  ;make the network look prettier
  repeat 1
  [
    layout-circle turtles (world-width / 2 - 1)
    ;layout-radial turtles links (turtle 0)
  ]
end


to sort-and-swap-turtles
  ; first sort the turtles by fear-score
  let x 0
  ; give turtles a label corresponding to their sorting order by fear score
  ;; and then randomly swap a small number of nodes
  foreach sort-on [fear-score] turtles [
    the-turtle -> ask the-turtle [
      set fear-score-sort x
      set x x + 1
    ]
  ]

   ; now iterate again and swap

   let i 0
   let my-position false
   let your-position false
   while [i <= count turtles][

    if (random-float 1.0 < percent-nodes-swapped)[
      ask one-of turtles [

        set my-position fear-score-sort
        ask one-of other turtles[
          set your-position fear-score-sort
          set fear-score-sort my-position
        ]

      set fear-score-sort your-position
      ]

    ]


   set i i + 1
   ]

end


to setup-homophily-small-world
  sort-and-swap-turtles ; sort turtles by fear-score with a bit of swapping
  set sort-idx-to-turtle-idx table:make
  set turtle-idx-to-sort-idx table:make
  layout-circle (sort-on [fear-score-sort] turtles) max-pxcor - 1
  let x 0
  ; give turtles a label corresponding to their sorting order by fear score
  ;; and then randomly swap a small number of nodes
  foreach sort-on [fear-score-sort] turtles [
    the-turtle -> ask the-turtle [
    table:put sort-idx-to-turtle-idx x who
    table:put turtle-idx-to-sort-idx who x

    set x x + 1
    ]
  ]


  ;Create a lattice
  ;; iterate over the turtles
  let n 0
  while [n < count turtles]
  [
    let i 1
    while [i <= int(average-node-degree / 2 )]
    [
      ;; make edges with the next (average-node-degree / 2) neighbors
      ;; this makes a lattice with average degree of average-node-degree
      make-edge turtle table:get sort-idx-to-turtle-idx n
      turtle (table:get sort-idx-to-turtle-idx ((n + i) mod count turtles))
      set i (i + 1)
    ]
    ;;50% chance of an extra link if average-node-degree mod 2 is non-zero (to compensate for odd average-node-degrees)
    if (random-float 1 < 0.5 and average-node-degree mod 2 > 0)
    [
      make-edge turtle n turtle ((n + int(average-node-degree / 2) + 1) mod count turtles)
    ]
    set n n + 1
  ]

 ;ask links [set rewired? false]
  ; in the spirit of a small world network, introduce some re-wiring
  repeat int(count links * percent-links-rewired)[

    let potential-edges links with [ not rewired? ]
    if any? potential-edges [
      ask one-of potential-edges [
        ;; "a" remains the same
        let node1 end1
        ;; if "a" is not connected to everybody
        if [ count link-neighbors ] of end1 < (count turtles - 1)
        [
          ;; find a node distinct from node1 and not already a neighbor of node1
          let node2 one-of turtles with [ (self != node1) and (not link-neighbor? node1) ]
          ;; wire the new edge
          ask node1 [ create-link-with node2 [ set color cyan  set rewired? true ] ]

          ;set number-rewired number-rewired + 1  ;; counter for number of rewirings

          ;; remove the old edge
          die
        ]
      ]
    ]
  ]

end



to setup-small-world-network

 layout-circle (sort turtles) max-pxcor - 1

  ;Create a lattice
  ;; iterate over the turtles
  let n 0
  while [n < count turtles]
  [
    let i 1
    while [i <= int(average-node-degree / 2 )]
    [
      ;; make edges with the next (average-node-degree / 2) neighbors
      ;; this makes a lattice with average degree of average-node-degree
      make-edge turtle n
              turtle ((n + i) mod count turtles)
      set i (i + 1)
    ]
    ;;50% chance of an extra link if average-node-degree mod 2 is non-zero (to compensate for odd average-node-degrees)
    if (random-float 1 < 0.5 and average-node-degree mod 2 > 0)
    [
      make-edge turtle n turtle ((n + int(average-node-degree / 2) + 1) mod count turtles)
    ]
    set n n + 1
  ]

  ;Rewire percent-links-rewired% of edges
  repeat int(count links * percent-links-rewired)[

    let potential-edges links with [ not rewired? ]
    if any? potential-edges [
      ask one-of potential-edges [
        ;; "a" remains the same
        let node1 end1
        ;; if "a" is not connected to everybody
        if [ count link-neighbors ] of end1 < (count turtles - 1)
        [
          ;; find a node distinct from node1 and not already a neighbor of node1
          let node2 one-of turtles with [ (self != node1) and (not link-neighbor? node1) ]
          ;; wire the new edge
          ask node1 [ create-link-with node2 [ set color cyan  set rewired? true ] ]

          ;set number-rewired number-rewired + 1  ;; counter for number of rewirings

          ;; remove the old edge
          die
        ]
      ]
    ]
  ]

end


;;;;;;;;;;;;;;;;
;; Clustering computations ;;
;;;;;;;;;;;;;;;;

to-report find-clustering-coefficient
  if (ticks = 0)[
  ; using the Watts-Strogatz definition
  ;; https://ccl.northwestern.edu/netlogo/docs/nw.html#nw:clustering-coefficient and https://arxiv.org/pdf/cond-mat/0303516.pdf p. 12
  set cc-nw mean [ nw:clustering-coefficient ] of turtles

]
  report cc-nw
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Path length computations ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report average-path-length

  if (ticks = 0 and behaviorspace-run-number = 0)[
    set apl-nw nw:mean-path-length
  ]
  if (behaviorspace-run-number != 0)[
    set apl-nw 0 ; do not calculate if we are in an experiment
  ]
  report apl-nw
end

;; connects the two turtles
to make-edge [node1 node2]
  ask node1 [ create-link-with node2  [
    set rewired? false
  ] ]
end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
1067
868
-1
-1
8.412
1
10
1
1
1
0
0
0
1
-50
50
-50
50
0
0
1
ticks
30.0

BUTTON
107
19
170
52
Go
Go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
28
19
92
52
Setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
0
117
198
150
Disease-Infection-Rate
Disease-Infection-Rate
0
1.0
0.1
0.01
1
NIL
HORIZONTAL

PLOT
1069
10
1363
218
Disease SIR Curves
Time
Agents
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"S" 1.0 0 -13345367 true "" ""
"I" 1.0 0 -2674135 true "" ""
"R" 1.0 0 -16777216 true "" ""

MONITOR
1127
222
1315
267
% Infected During Simulation
percent-infected
3
1
11

SLIDER
3
159
175
192
Infection-Duration
Infection-Duration
0
20
6.0
1
1
NIL
HORIZONTAL

TEXTBOX
8
98
158
116
Pathogen Parameters:
11
0.0
1

TEXTBOX
1387
22
1547
103
Color Code\n\nBlue - Susceptible\nRed - New Adopter\nWhite - Old Adopter
12
0.0
1

SLIDER
2
771
174
804
average-node-degree
average-node-degree
0
15
6.0
1
1
NIL
HORIZONTAL

SLIDER
0
212
172
245
num-turtles
num-turtles
0
5000
1000.0
1
1
NIL
HORIZONTAL

CHOOSER
2
813
194
858
Network-Type
Network-Type
"Spatially Clustered" "Random" "Small World" "By Fear-Score" "Small World Fear-Score"
0

SLIDER
3
859
184
892
percent-links-rewired
percent-links-rewired
0
1
0.2
0.01
1
NIL
HORIZONTAL

SWITCH
7
278
121
311
Fear-On?
Fear-On?
0
1
-1000

TEXTBOX
2
747
152
765
Network parameters:
11
0.0
1

TEXTBOX
10
259
160
277
Fear Parameters:
11
0.0
1

SLIDER
2
361
176
394
Fear-Infection-Rate
Fear-Infection-Rate
0
1
0.05
0.01
1
NIL
HORIZONTAL

SLIDER
2
395
174
428
hide-duration
hide-duration
0
14
4.0
1
1
NIL
HORIZONTAL

BUTTON
107
53
201
86
Go (Once)
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
7
311
232
356
Fear-Mechanism
Fear-Mechanism
"Complex Contagion" "FTA" "Complex Contagion-Gamma" "FTA-Mean" "FTA-v2"
2

SLIDER
3
488
175
521
fear-threshold
fear-threshold
1
10
1.0
1
1
NIL
HORIZONTAL

TEXTBOX
5
447
176
475
# of Fear Transmissions to Hide (set to 1 for simple contagion):
11
0.0
1

TEXTBOX
0
609
176
637
For Fear-Mechanism = \"FTA\"
11
0.0
1

SLIDER
2
670
174
703
satisficing-radius
satisficing-radius
0
0.5
0.01
0.01
1
NIL
HORIZONTAL

SWITCH
0
632
126
665
satisficing?
satisficing?
1
1
-1000

PLOT
1390
436
1590
586
Fear-Score Distribution
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"set-plot-pen-mode 1\nset-plot-pen-interval 0.025" ""
PENS
"pen-0" 1.0 0 -7500403 true "" "histogram [fear-score] of turtles"

PLOT
1114
439
1314
589
Number Hiding
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"H" 1.0 0 -16777216 true "" "plot count turtles with [hiding? = True]"

SLIDER
1511
715
1683
748
fear-lambda
fear-lambda
0.1
1
1.0
0.1
1
NIL
HORIZONTAL

MONITOR
1592
436
1705
481
Mean Fear Score
mean [fear-score] of turtles
5
1
11

TEXTBOX
1396
677
1546
705
Initial Fear Score Configuration for \"FTA\"
11
0.0
1

MONITOR
1134
288
1272
333
Day w/ Max Infected
day-of-max-infected
17
1
11

PLOT
1128
629
1328
779
Count Fearful
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count turtles with [fear-state = \"I\"]"

INPUTBOX
1144
837
1293
897
random_seed_num
1.0
1
0
Number

SWITCH
1314
865
1463
898
Random-Seed?
Random-Seed?
0
1
-1000

SWITCH
1368
732
1518
765
color-patches?
color-patches?
1
1
-1000

SLIDER
3
707
175
740
step-increment
step-increment
0
0.2
0.005
0.0025
1
NIL
HORIZONTAL

SLIDER
1593
43
1765
76
init-infected
init-infected
1
10
1.0
1
1
NIL
HORIZONTAL

SLIDER
1541
766
1713
799
beta-dist-denom
beta-dist-denom
0
10
1.0
1
1
NIL
HORIZONTAL

SLIDER
164
273
336
306
restore-hiding
restore-hiding
0
500
0.0
5
1
NIL
HORIZONTAL

MONITOR
1592
482
1694
527
Var Fear Score
variance [fear-score] of turtles
3
1
11

CHOOSER
1
536
193
581
type-of-fear-distribution
type-of-fear-distribution
"Gamma" "Anchor" "Custom Normal" "Custom Normal v2" "Uniform" "Truncated Normal"
0

MONITOR
1458
311
1631
356
Average Path Length
average-path-length
2
1
11

MONITOR
1456
264
1634
309
Clustering Coefficient
find-clustering-coefficient
2
1
11

SLIDER
223
890
429
923
percent-nodes-swapped
percent-nodes-swapped
0
1
1.0
0.01
1
NIL
HORIZONTAL

PLOT
1268
282
1468
432
Avg Fear Score of Infected
NIL
NIL
0.0
10.0
0.0
0.25
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

CHOOSER
542
910
689
955
Release-Mechanism
Release-Mechanism
"fixed-duration" "probabilistic"
0

SLIDER
1545
830
1740
863
prob-release-damping
prob-release-damping
0
200
150.0
1
1
NIL
HORIZONTAL

SWITCH
1594
628
1755
661
export-network?
export-network?
1
1
-1000

SWITCH
114
236
411
269
agents-can-spread-fear-while-hiding?
agents-can-spread-fear-while-hiding?
0
1
-1000

@#$#@#$#@
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment-01-fear-increasing-spatially-clustered" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Spatially Clustered&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.08"/>
      <value value="0.09"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-02-fear-decreasing-spatially-clustered" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Spatially Clustered&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.08"/>
      <value value="0.09"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-03-fear-increasing-small-world" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-04-fear-decreasing-small-world" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-05-fear-increasing-small-world-vary-clustering-and-swapping" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>percent-infected</metric>
    <metric>fear-score-of-infected</metric>
    <metric>fear-score-of-infected-cumulative</metric>
    <metric>count-hiding</metric>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World Fear-Score&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="percent-links-rewired" first="0.1" step="0.01" last="0.2"/>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;FTA&quot;"/>
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="percent-nodes-swapped" first="0" step="0.25" last="1"/>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-06-fear-decreasing-small-world-vary-clustering-and-swapping" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>percent-infected</metric>
    <metric>fear-score-of-infected</metric>
    <metric>fear-score-of-infected-cumulative</metric>
    <metric>count-hiding</metric>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World Fear-Score&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="percent-links-rewired" first="0.1" step="0.01" last="0.2"/>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;FTA&quot;"/>
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="percent-nodes-swapped" first="0" step="0.25" last="1"/>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-07-fear-decreasing-small-world-probabilistic-release" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;probabilistic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="150"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="calibrate-prob-release-damping" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="150"/>
    <metric>mean avg-hide-time</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;probabilistic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <steppedValueSet variable="prob-release-damping" first="100" step="10" last="200"/>
  </experiment>
  <experiment name="experiment-01a-fear-increasing-spatially-clustered-complex-contagion" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Spatially Clustered&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.08"/>
      <value value="0.09"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-03a-fear-increasing-small-world-complex-contagion" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-02a-fear-decreasing-spatially-clustered-complex-contagion" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Spatially Clustered&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.08"/>
      <value value="0.09"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-04a-fear-decreasing-small-world-complex-contagion" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-08-one-shot-example-Contagion" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="150"/>
    <metric>[(list who ticks disease-state fear-state fear-score tick-infected dend)] of turtles</metric>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World Fear-Score&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-08-one-shot-example-FTA" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-world "../test/one-shot-FTA.csv"</final>
    <timeLimit steps="150"/>
    <metric>[(list who ticks disease-state fear-state fear-score tick-infected dend)] of turtles</metric>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World Fear-Score&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-09-fear-decreasing-small-world-limited-uptake-under-complex-contagion" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-00-show-sample-network" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>export-world "../data/core/sample-network.csv"</final>
    <timeLimit steps="40"/>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Spatially Clustered&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-network?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-03-fear-increasing-small-world-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-04-fear-decreasing-small-world-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-02-fear-decreasing-spatially-clustered-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Spatially Clustered&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.08"/>
      <value value="0.09"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-01-fear-increasing-spatially-clustered-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Spatially Clustered&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.08"/>
      <value value="0.09"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-05-fear-increasing-small-world-vary-clustering-and-swapping-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>percent-infected</metric>
    <metric>fear-score-of-infected</metric>
    <metric>fear-score-of-infected-cumulative</metric>
    <metric>count-hiding</metric>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World Fear-Score&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="percent-links-rewired" first="0.1" step="0.01" last="0.2"/>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;FTA&quot;"/>
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="percent-nodes-swapped" first="0" step="0.25" last="1"/>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-06-fear-decreasing-small-world-vary-clustering-and-swapping-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>percent-infected</metric>
    <metric>fear-score-of-infected</metric>
    <metric>fear-score-of-infected-cumulative</metric>
    <metric>count-hiding</metric>
    <enumeratedValueSet variable="random_seed_num">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World Fear-Score&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="percent-links-rewired" first="0.1" step="0.01" last="0.2"/>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;FTA&quot;"/>
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="percent-nodes-swapped" first="0" step="0.25" last="1"/>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-01a-fear-increasing-spatially-clustered-complex-contagion-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Spatially Clustered&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.08"/>
      <value value="0.09"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-02a-fear-decreasing-spatially-clustered-complex-contagion-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Spatially Clustered&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.08"/>
      <value value="0.09"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-03a-fear-increasing-small-world-complex-contagion-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-04a-fear-decreasing-small-world-complex-contagion-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;fixed-duration&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-07-fear-decreasing-small-world-probabilistic-release-sensitivity" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-disease-I</metric>
    <metric>count-fearful</metric>
    <metric>count-hiding</metric>
    <metric>percent-infected</metric>
    <enumeratedValueSet variable="beta-dist-denom">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing-radius">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-lambda">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-infected">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fear-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="step-increment">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-duration">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Type">
      <value value="&quot;Small World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="satisficing?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-fear-distribution">
      <value value="&quot;Anchor&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-nodes-swapped">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-turtles">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Infection-Duration">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Random-Seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Disease-Infection-Rate">
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-patches?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-links-rewired">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Mechanism">
      <value value="&quot;Complex Contagion-Gamma&quot;"/>
      <value value="&quot;FTA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-Infection-Rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fear-On?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_seed_num">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="restore-hiding">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Release-Mechanism">
      <value value="&quot;probabilistic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-release-damping">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agents-can-spread-fear-while-hiding?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
