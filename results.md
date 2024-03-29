# Ratio calculations using search.py (?)

**div** -- number of divisions on each side, i.e. genus = div^dim

**cdist** -- maximum allowed l1-distance between cubes, i.e., codimension of intersection (* = no restrictions)

# 2D curves

## 2D monofractals

|div|entr-exit      |cdist|#upaths|#paths|l1   |l2   |linf |
|---|---------------|-----|-------|------|-----|-----|-----|
|2  |side           |*    |1      |      |**9.000**|6.000|6.000|
|3  |side           |*    |5      |      |9.800|5.666|**4.500**|
|3  |diag           |*    |2      |      |10.00|5.666|5.333|
|4  |side           |1    |5      |      |9.000|5.909|5.4545|
|4  |side           |*    |162    |298   |9.000|5.909|5.294|
|5  |side           |1    |43     |      |10.00|**5.589**|5.333|
|5  |side           |*    |24850  |49700 |5.589|     |     |
|5  |diag           |1    |11     |      |10.66|6.0  |5.77 |
|5  |diag           |*    |659    |2592  |6.0  |     |     |
|5  |(0,0)->(1/2,1) |*    |557    |      |11.85|6.535|6.25 |
|6  |side           |1    |897    |      |9.000|5.62 |5.294|
|6  |side           |*    |       |      |     |     |     |
|6  |(0,0)->(1/2,1) |1    |1476   |      |     |6.5  |     |

## 2D bifractal curves, div=2

cdist=* for all rows.

|gate1           |gate2            |#paths|l1  |l2  |linf|
|----------------|-----------------|------|----|----|----|
|(0,0)->(0,1/2)  |(0,0)-->(1/2,1)  |8     |13.3|6.66|5.0 |
|(0,0)->(0,1)    |(0,0)->(0,1)     |1     |9.0 |6.0 |6.0 |
|(0,0)->(0,1)    |(0,0)-->(1/2,1)  |3     |12.0|6.0 |5.0 |
|(0,0)->(0,1)    |(0,0)-->(2/3,1)  |2     |12.0|6.0 |5.66|
|(0,0)->(0,1)    |(0,0)-->(1,1)    |5     |12.0|6.0 |5.14|
|(0,0)->(0,1)    |(0,1/2)-->(0,1/2)|2     |9.0 |6.0 |6.0 |
|(0,0)->(0,1)    |(0,1/2)-->(1/2,0)|2     |9.0 |6.0 |6.0 |
|(0,0)->(0,1)    |(0,1/2)-->(1,1/2)|2     |9.0 |6.0 |6.0 |
|(0,1/3)->(1/3,1)|(0,1/3)-->(1,1/3)|1     |9.0 |**5.0**|5.0 |

## 2D bifractal curves, div=3

|gate1           |gate2            |#paths|l1  |l2  |linf|
|----------------|-----------------|------|----|----|----|
|(0,0)->(0,1/3)  |(0,0)->(0,1)     |90    |11.3|6.25|5.0 |
|(0,0)->(0,1/2)  |(0,0)->(0,1)     |81    |11.3|6.25|5.0 |
|(0,0)->(0,1)    |(0,0)->(0,1)     |15    |    |6.25|9.8 |
|(0,0)->(0,1)    |(0,0)->(1/2,1)   |864   |11.3|6.25|5.0 |
|(0,0)->(0,1)    |(0,0)->(2/3,1)   |2744  |9.8 |5.66|4.5 |
|(0,0)->(0,1)    |(0,1/3)->(0,2/3) |45    |9.8 |5.66|4.5 |
|(0,0)->(0,1)    |(0,1/3)->(1/3,1) |75    |9.8 |5.66|4.5 |
|(0,0)->(0,1)    |(0,1/3)->(1,1/3) |25    |9.8 |5.66|4.5 |
|(0,0)->(0,1)    |(0,0)->(1,1)     |23010 |    |?   |    |
|(0,0)->(1/4,1)  |(0,0)->(3/4,1)   |2560  |9.8 |5.8 |5.33|
|(0,0)->(1,1)    |(0,0)->(1,1)     |3     |    |?   |    |
|(0,1/4)->(1/4,1)|(0,1/4)->(1,1/4) |2     |10.6|5.66|4.5 |
|(0,1/4)->(1/4,1)|(0,1/4)->(1,3/4) |4     |9.8 |5.0 |4.5 |

# 3D curves

## 3D monofractal curves, div=2

|entr-exit              |cdist|#paths|l1   |l2   |linf |
|-----------------------|----|------|-----|-----|-----|
|(0,0,0)->(0,0,1)       |*   |29    |**89.744**|22.86|**12.4**|
|(0,0,0)-->(0,1,1)      |*   |149   |99.13|21   |13   |
|(0,0,0)-->(0,1/2,1)    |*   |2758  |102.9|25.6 |16.59|
|(0,0,0)->(1/2,1/2,1)   |*   |4     |180.45|40.2 |28.0|
|(0,0,1/3)-->(1/3,1,1)  |*   |16    |98.0 |23.6 |16.59|
|(0,1/3,1/3)-->(1/3,1/3,1)|* |1     |89.76|**18.56**|14.0 |

## 3D bifractal curves, div=2

Gates: 401 total;

* 0/4 vertex: 202
(add means additional gate entrance/exit coordinates in {0,1})
  * 0/4 add: 35, with l2 below 20: 19; best l2: **16.991**; gates: "(0,1/3,1/3)->(1/3,0,2/3)", "(0,1/3,1/3)->(1/3,1/3,1)"
  * 1/4 add: 2, with l2 below 20: 0
  * 2/4 add: 47, with l2 below 20: 0
  * 3/4 add: 58
  * 4/4 add: 60

* 1/4 vertex: 135
* 2/4 vertex: 31
* 3/4 vertex: 28
* 4/4 vertex: 5
