```{include} ./header.md
```

# Component map

```{tableofcontents}
```

Meshing scripts use base components. The list of predefined elements is shown with `d_rail('nmap')`

## Data-base architecture

Current tests done  in {m}`t_exp19('script')`


List of contents  <a href="matlab:d_rail('meshdb;')">{m}`d_rail('meshdb')`</a>.

Current tests <a href="matlab:sdtweb(t_exp19,'ref')">{m}`t_exp19('ref')`</a>.


:::{dropdown} <a href="matlab:d_rail('meshdb;')">![](./_icons/run16.png){m}`d_rail('meshdb')`</a>. Expand to display source code.
```matlab
 d_rail('meshdb')
```
:::

:::{list-table} Initial database in Map:Sections
:header-rows: 0
:widths: auto
*   - ![](./_images/U30_e.png)
    - ![](./_images/U30_ra.png)
    - ![](./_images/U30_ra+be.png)
:::
<img src="./_images/U30_ra+e.png" alt="U30_ra+e" width="200"/>
<img src="./_images/U30_ra+s.png" alt="U30_ra+s" width="200"/>
<img src="./_images/U30_ra+sbe.png" alt="U30_ra+sbe" width="200"/>
<img src="./_images/U30_ra+se.png" alt="U30_ra+se" width="200"/>
<img src="./_images/Wheel.png" alt="Wheel" width="200"/>
<img src="./_images/WheelCut.png" alt="WheelCut" width="200"/>

- Map:Sections contains 
  - U30.ra : 2D section of U30 rail 
  - U30.ra+s : rail + sleeper 2D section 
  - Ec41 : list of positions for a Bracket   
  - xxx 


:::{list-table} Initial database in Map:Sections
:widths: auto
:header-rows: 1
*   - key
    - ToolTip
    - content
*   - Ec0
    - 
    - {'ra',-300,[],15;
'ra+s',-127.5,[],15;
'ra',127.5,[],15;
'',300,[],15}
*   - Ec41
    - 4 bolts, also called E120
    - struct('ToolTip','4 bolts, also called E120', ...
'li',{{'ra+e',-289,[],15;
'ra+se',-274,[],15;
'ra+sbe',-259,[],15;
'ra+se',-211,[],__}})
*   - Ec42
    - 
    - {'ra+e',-289,[],15;
'ra+e',-274,[],15;
'ra+be',-259,[],15;
'ra+e',-211,[],15,__}
*   - Ec61
    - 6 bolts
    - struct('ToolTip','6 bolts', ...
'li',{{'ra+e',-459,[],15;
'ra+e',-444,[],15;
'ra+be',-429,[],15;
'ra+e',-381,[],15,__}})
*   - Ec62
    - 
    - {'ra+e',-459,[],15;
'ra+e',-444,[],15;
'ra+sbe',-429,[],15;
'ra+se',-381,[],15,__}
*   - U30.e
    - 
    - 116 Node,81 Elt
*   - U30.ra
    - 
    - 152 Node,109 Elt
*   - U30.ra+be
    - 
    - 2870 Node,1483 Elt
*   - U30.ra+e
    - 
    - 294 Node,190 Elt
*   - U30.ra+s
    - 
    - 236 Node,184 Elt
*   - U30.ra+sbe
    - 
    - 3290 Node,1783 Elt
*   - U30.ra+se
    - 
    - 378 Node,265 Elt
*   - Wheel
    - 
    - 14560 Node,10488 Elt
:::


## Rails

- U30 is the rail type (2D section)
- bU30 xxx beam, see dynavoie 
- cu : cube representation of rail head, ScldCS xxxeb


## Sleepers

- xxx


## Brackets

d_rail('nmap.Map:Bracket')


## Pads 

- xxx


## Mat : named materials


d_rail('nmap.Map:Mat')



## Wheels


