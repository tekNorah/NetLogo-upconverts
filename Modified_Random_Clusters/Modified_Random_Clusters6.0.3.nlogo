;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Modified random clusters method for landscape pattern simulation - NetLogo  ;;
;;                                                                              ;;
;;  Code licenced by James D.A. Millington (http://www.landscapemodelling.net)  ;;
;;  under a Creative Commons Attribution-Noncommercial-Share Alike 3.0          ;;
;;  Unported License (see http://creativecommons.org/licenses/by-nc-sa/3.0/)    ;;
;;                                                                              ;;
;;  This code is based on the method by Saura & Martinez-Millan (2000)          ;;
;;  Landscape Ecology 15 661-678 http://dx.doi.org/10.1023/A:1008107902848      ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



patches-own
[
  cover        ;; cover of patch
  cluster_ID   ;; id of cluster patch belongs to
  cluster      ;; holds a patch which is that cluster's "leader" - arbitrarily chosen
]


globals
[
  cluster-area-list  ;for plotting purposes
]




;---------------
to generate-landscape

  set cluster-area-list []

  set_random_seeds                                 ;; randomly distribute cluster seeds (Step 1 in Saura & Martinez-Millan 2000 p.664)

  set-cover                                        ;; identify and label clusters (Steps 2 & 3 in above)
  fill_landscape                                   ;; fill remaining patches to dominant neighbour (Step 4 in above)

  ask patches [ set pcolor (cover * 10) + 5 ]      ;; color patches according to cover types

  plotting                                         ;; plot bar plots

end
;---------------




;---------------
to set-cover
  loop
  [
    let seed one-of patches with [ ( cluster = nobody and pcolor = blue ) ]   ;; pick a random seed we haven't labelled yet
    if seed = nobody                                                          ;; if we have labelled all seeds stop
      [ stop ]

    ask seed
    [
      set cluster self                     ;; make the seed the "leader" of a new cluster by assigning itself to its own cluster
      set cover ( random Number-of-Types )       ;; assign this "leader" a random cover
      grow-cover_cluster                   ;; then grow-cover_cluster to find the rest of the cluster

      set cluster-area-list fput count patches with [cluster = myself] cluster-area-list
    ]

  ]
end
;---------------




;---------------
to grow-cover_cluster

  without-interruption              ;; we use without-interruption here so each patch only gets added to the cluster once
  [
    let neighbours nobody

    ifelse(Neighbourhood = "Moore")
    [ set neighbours neighbors with [ cluster = nobody and pcolor = [pcolor] of myself ] ]
    [ set neighbours neighbors4 with [ cluster = nobody and pcolor = [pcolor] of myself ] ]

    ask neighbours
    [
      set cluster [cluster] of myself    ;;assign neighbouring patch to seed patch cluster
      set cover [cover] of myself        ;;assign neighbouring patch to seed patch cover
      set cluster_ID [cluster_ID] of myself
      grow-cover_cluster                 ;;recursive call!
    ]
   ]

end
;---------------





;---------------
to set_random_seeds

  let cluster-counter 0

  ask patches
  [
    set cover nobody
    set cluster_ID nobody
    set cluster nobody
  ]

  ask patches
  [
    ifelse ( (random-float 1) < cluster-probability )  ;set patch as seed if rand no is less than cluster_probability
    [
      set pcolor blue                                ;seeds are blue
      set cluster_ID cluster-counter
      set cluster-counter cluster-counter + 1
    ]
    [ set pcolor red ]                               ;non-seeds are red
  ]

end
;---------------



;---------------
to fill_landscape                  ;; assign patches not initially assigned (to dominant cover in neighbourhood)

  loop
  [
    if ( not any? patches with [ cluster = nobody ] )     ;;if all patches have been assigned to an cover-cluster
    [ stop ]

    ask one-of patches with [ cluster = nobody ]  ;;pick a patch that has not been assigned to an cover-cluster yet
    [
      let neighbours nobody

      ifelse(Neighbourhood = "Moore")
      [ set neighbours neighbors ]
      [ set neighbours neighbors4 ]


      ifelse ( any? neighbours with [ cluster != nobody ] )  ;; check if there are any assigned patches in neighbourhood
      [

        let covers []

        ask neighbours with [ cluster != nobody ]
        [
          set covers fput cover covers    ;;ask neighbours to add their covers to the list
        ]

        let unique-covers remove-duplicates covers    ;;create a list of unique covers

        let max-cover-count -1                 ;the number of neighbours with the maximum cover
        let max-cover -1                       ;the maximum cover


        ifelse(length unique-covers > 1)
        [
          ;if there is more than one unique-cover
          foreach unique-covers                  ;for each of the unique covers
          [ [?1] ->
            let occ occurrences ?1 covers          ;count how many neighbours had this cover

            ifelse(occ > max-cover-count)        ;if the count is greater than the current maximum count
            [
              set max-cover ?1                    ;set this as the dominant cover
              set max-cover-count occ            ;update the current maximum count
            ]
            [
              if(occ = max-cover-count)          ;otherwise, if the count is equal to the current maximum count
              [
                let rand random-float 1
                if(rand < 0.5) [ set max-cover ?1 ]  ;randomly pick one of the two covers to be the new dominant cover
              ]
            ]
          ]
        ]

        [
          ;otherwise just set the max-cover to the only unique-cover
          set max-cover first unique-covers
        ]


        let p one-of neighbours with [ cover = max-cover ]    ;;assign qualities of one of dominant neighbours to patch (no wrap!)
        set pcolor [pcolor] of p
        set cluster [cluster] of p
        set cover [cover] of p
        set cluster_id [cluster_id] of p
      ]
      [                                                               ;; if no assigned agents in neighbourhood assign an agent at random
        set cover ( random Number-of-Types)
        set cluster self                                              ;; remember to start a new cluster!
        set cluster_ID ([cluster_ID] of max-one-of patches [cluster_ID]) + 1
      ]
    ]
  ]

end
;---------------





;; count the number of occurrences of an item in a list
;---------------
to-report occurrences [x the-list]
  report reduce
    [ [?1 ?2] -> ifelse-value (?2 = x) [?1 + 1] [?1] ] (fput 0 the-list)
end
;---------------



;---------------
to plotting

  set-current-plot "Cluster-Distribution"
  clear-plot
  histogram cluster-area-list

  set-current-plot "Cover-By-Area"
  clear-plot
  let cover-counts []
  let c 0
  while[c < Number-of-Types]
  [
    set cover-counts fput 0 cover-counts
    set c c + 1
  ]
  ask patches [ set cover-counts replace-item cover cover-counts ((item cover cover-counts) + 1) ]

  set c 0
  while[c < Number-of-Types]
  [
    plotxy c item c cover-counts
    set c c + 1
  ]



  set-current-plot "Count-Clusters-By-Cover"
  clear-plot

  set cover-counts []
  set c 0
  while[c < Number-of-Types]
  [
    set cover-counts fput 0 cover-counts
    set c c + 1
  ]

  let max-cluster-ID 0
  ask max-one-of patches [cluster_ID][ set max-cluster-ID cluster_ID ]

  let id 0
  while[id <= max-cluster-ID]
  [
    let p one-of patches with [cluster_ID = id]

    if(p != nobody)
    [
      ask p
      [
        set cover-counts replace-item cover cover-counts ((item cover cover-counts) + 1)
      ]
    ]
    set id id + 1
  ]

  print cover-counts

  set c 0
  while[c < Number-of-Types]
  [
    plotxy c item c cover-counts
    set c c + 1
  ]

end
;---------------



;---------------
to export-map

  file-close-all
  if(file-exists? "export.asc") [file-delete "export.asc"]
  file-open "export.asc"


  file-type "ncols   "
  file-print world-width
  file-type "nrows   "
  file-print world-height
  file-print "xllcorner  0"
  file-print "yllcorner  0"
  file-print "cellsize   1"
  file-print "NODATA_value  -9999"

  let x min-pxcor
  let y max-pycor

  while[y >= min-pycor]
  [
    while[x <= max-pxcor]
    [
      ifelse(x = max-pxcor)
      [
        ask patches with [pxcor = x and pycor = y] [ file-print cover ]
      ]
      [
        ask patches with [pxcor = x and pycor = y] [ file-type cover file-type " "]
      ]

      set x x + 1
    ]

    set x min-pxcor
    set y y - 1
  ]

  file-close

end
;---------------
@#$#@#$#@
GRAPHICS-WINDOW
221
66
637
483
-1
-1
8.0
1
10
1
1
1
0
1
1
1
-25
25
-25
25
0
0
1
ticks
30.0

SLIDER
28
109
200
142
number-of-types
number-of-types
1
10
5.0
1
1
NIL
HORIZONTAL

SLIDER
27
198
199
231
cluster-probability
cluster-probability
0.05
1
0.5
0.05
1
NIL
HORIZONTAL

BUTTON
28
70
190
103
NIL
generate-landscape
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
25
301
164
346
Neighbourhood
Neighbourhood
"von Neumann" "Moore"
0

PLOT
648
61
848
211
Cluster-Distribution
Area
Count
1.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

PLOT
648
215
848
365
Count-Clusters-by-Cover
Cover
Count
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

PLOT
649
368
849
518
Cover-By-Area
Cover
Area
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

BUTTON
27
420
135
453
NIL
export-map
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
31
145
181
177
Number of types to generate in the landscape
13
0.0
1

TEXTBOX
31
233
208
297
Probability for percolation map generation (see Step A Saura & Martinez-Millan 2000)
13
0.0
1

TEXTBOX
28
355
204
419
Neighbourhood rule in percolation (von Neumann is neighbors4)
13
0.0
1

TEXTBOX
29
464
198
512
ESRI ascii grid format (to \\\"export.asc\\\", overwritten if exists; does not work online)
13
0.0
1

TEXTBOX
223
10
660
50
Modified random clusters method for landscape pattern simulation
15
0.0
1

TEXTBOX
231
32
658
63
NetLogo code by  by James D.A. Millington (http://www.landscapemodelling.net) based on method by Saura & Martinez-Millan (2000, httpp://dx.doi.org/10.1023/A:1008107902848)
10
0.0
1

@#$#@#$#@
## WHAT IS IT?

A implementation of the modified random clusters method for landscape pattern simulation described by Saura and Martinez-Millan (2000) and implemented in SIMMAP (http://www2.montes.upm.es/personales/saura/software.html#simmap). Used previously in Millington et al. (2008, http://jasss.soc.surrey.ac.uk/11/4/4.html)

## HOW IT WORKS

Follows the Steps A-D described in Saura and Martinez-Millan (2000, p.664-665)

## HOW TO USE IT

Select the number of land cover (or habitat, etc.) types to be generated in the map. Select the cluster-probability to determine amount of clustering. Select whether von Neumann or Moore neighbourhood rule should be used. Click 'generate-landscape' button. Once landscape pattern is generated the map can be exported to ESRI ascii grid format by clicking "export-map".

## THINGS TO TRY

Examine how the size and number of patches produced varies for different values of number-of-types, cluster-probability and neighbourhood type.

## EXTENDING THE MODEL

Feel free to use the code in your own models.

## RELATED MODELS

Similar code used in Millington, J.D.A., Romero-Calcerrada, R., Wainwright, J. and Perry, G.L.W. (2008) An agent-based model of Mediterranean agricultural land-use/cover change for examining wildfire risk Journal of Artificial Societies and Social Simulation 11(4) 4 http://jasss.soc.surrey.ac.uk/11/4/4.html

## CREDITS AND REFERENCES

Modified random clusters method described in: Landscape patterns simulation with a modified random clusters method, Santiago Saura and Javier Mart�nez-Mill�n (2000) Landscape Ecology 15 661-678 http://dx.doi.org/10.1023/A:1008107902848
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
NetLogo 6.0.3
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
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
