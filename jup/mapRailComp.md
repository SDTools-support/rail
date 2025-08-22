# Component map

Meshing scripts use base components. The list of predefined elements is shown with `d_rail('nmap')`

 - [Map:Sections](d_rail.sections) contains mesh components that are a combined to form the track
 - [Map:Slice](d_rail.slice) correspond to 

Current tests done in {m}`t_exp19('script')`


List of contents  <a href="matlab:d_rail('meshdb;')">{m}`d_rail('meshdb')`</a>.

Current tests <a href="matlab:sdtweb(t_exp19,'ref')">{m}`t_exp19('ref')`</a>.


<a href="matlab:d_rail('meshdb;')"> .</a>

````{dropdown} <a href="matlab:gartid('tutoBy-s2;')">![](../_images/run16.png)Run </a>. Expand to display source code.
```matlab
 d_rail('meshdb')
```
````


(d_rail.sections)=
## Map:Sections map of mesh components

Mesh building is performed by combining a number of mesh components which should be loaded into the project. xxx

````{list-table} Armament labels 
:widths: 40 20 20 20
:header-rows: 0

* - Rail :2D alone `ra`, with sleeper `ra+s` (optional)
  - ![](../_images/U30_ra.png)
  - ![](../_images/U30_ra+s.png)
  - 
* - Fishplate no bolt (*éclisse* `e` 2D)
  - ![](../_images/U30_e.png)
  - ![](../_images/U30_ra+e.png)
  - ![](../_images/U30_ra+se.png)

* - Fishplate bolt (RD)
  - ![](../_images/U30_ra+be.png)
  - ![](../_images/U30_ra+sbe.png)
  -
* - Wheel
  - ![](../_images/Wheel.png)
  - ![](../_images/WheelCut.png)
  -
````

Profiles for [wheel](cntc.prw) and [rail](cntc.prr) profiles can be xxx

(d_rail.slices)=
## Map:Slices map of pre-mesh slice parameters


````{list-table} Initial database in Map:Slices
:widths: 10 20 70 
:header-rows: 1
*   - key
    - ToolTip
    - content
*   - Ec0
    - Rail (no bolt)
    - `{'ra',-300,[],15;'ra+s',-127.5,[],15;'ra',127.5,[],15;'',300,[],15}`
*   - Ec41
    - 4 bolts, also called E120
    - `struct('ToolTip','4 bolts, also called E120', 'li',{{'ra+e',-289,[],15;'ra+se',-274,[],15;'ra+sbe',-259,[],15;'ra+se',-211,[],xxx}})`
````


(Rails)=
## Rails

- U30 is the rail type (2D section)
- bU30 xxx beam, see dynavoie 
- cu : cube representation of rail head, ScldCS xxxeb


(Sleepers)=
## Sleepers

- xxx


(Brackets)=
## Brackets

d_rail('nmap.Map:Bracket')


## Pads 

- xxx


## Mat : named materials


d_rail('nmap.Map:Mat')



## Wheels


