# urn for rail wheel definition

```{tableofcontents}
```

## Simple example

Meshing scripts use base components. The list of predefined elements is show with `d_rail('nmap')` 
 
|key      |def     |ToolTip                           |fmt|
|---      |---     |---                               |---|
|gr.All   |        |Experiment parameters             |   |
|gr.s1    |        |Slice definition                  |   |
|rail     |        |Rail Urn, eg U60                  |%s |
|pad      |        |Pad Urn, eg PadFuSn{io4}          |%s |
|sleeper  |        |Sleeper Urn, eg mass{4}           |%s |
|sub      |        |Sub-structure Urn, eg spring{k,c} |%s |
|gr.Track |        |Track definition parameters       |   |
|top      |-700 300|contact refine start end          |%g |
|topCoarse|        |do not refine coarse              |31 |
|sw       |.6      |slice width                       |%s |
|quad     |        |use quadratic elements            |31 |
|Lc       |50      |characteristic length             |%g |
|gr.Wheel |        |Wheel refinement parameters       |   |
|Wref     |coarse1 |wheel refinement strategy         |%s |
|gr.Ctc   |        |Contact parameters                |   |
|Kc       |-23e9   |contact stiffness (>0 lin,<0 sqrt)|%ug|
|delta    |20      |contact lowpass                   |%ug|

The mesh is obtained using by extruding and combining base elements separated by `:`. For example 
 `U30{5,Gc,tc-350_400}:Ec41{Air}:W2{XaZa}` (from `rail19('nmap.Trk21ref.Ref21')`) combines

- `U30{5,Gc,tc-350_400}`  corresponds to the base rail nomenclature where 
  -  `U30` is the [rail](mapRailComp.md#Rails) type 
  -  `Ns=5` is the number of sleepers before and after the gap,  xxxdynavoie {n\_slice} 
  - {s}`Gc` specificity the gradient type {m}`TGrag=c`. See the tag {m}`rail19('Mesh.TGrad')` for implementation. 
  -  {s}`tc-350\_400` t is for top, tc indicates RM.topCoarse=1,  RM.top limits are [-350 500], see \ser{refinerail}
 - `Ec41{Air}` specifies a bracket (*éclisse*) configuration of type `Ec41` with an Air gap. Interpretation of associated configurations in done in {m}`rail19('MailRailDb')`. A list of configurations is given in the {m}`Ec` entries of {m}`rail19('nmap')`. See for example {m}`rail19('nmap.Ec41')` the 4 bolt configuration. See figure \ref{diffeclisse11} for supported list. 
- `Wa{XaZu}` specifies the wheel type and trajectory nature in horizontal x and vertical z directions. {m}`xa` stands for acceleration trajectory obtained by specifying a distributed load, {m}`xu` stands for an enforced displacement strategy. Possibly use {m}`W0` for a point load, see \ser{WheelTraj}


## URN interpretation

Interpretation of the uniform name is performed in `d_rail('nameToMeshRO')`.

When generating the mesh a number of sections are considered (building of sections is done in {m}`rail19('Build\_sections')`) and show in figure~\ref{sections}.
%
\begin{itemize}
\item {s}`ra+s` contains a rail + sleeper cut (thus including the pad mesh)
\item {s}`ra` rail only 
\item {s}`ra+e` rail + bracket ({\em éclisse}), {m}`se` sleeper and bracket, {m}`sbe` sleeper, bolt, bracket, {m}`rail+e:Air` rail gap (using Air as material) + bracket.  
\end{itemize}

\begin{figure}[H]
		\centering
    	\includegraphics[width=.49\textwidth]{modelschema3}
       \includegraphics[width=.49\textwidth]{boulon}
    \caption{Named sections {m}`ra` rail, {m}`e` bracket, {m}`ra+e` ,{m}`ra+se`, {m}`ra+sbe` bolt over a sleeper, {m}`ra+be` bolt not over a sleeper\label{sections}}
\end{figure}

The simulation configuration nomenclature is currently not used. 



The material nomenclature can be found using {m}`d_rail('nmap.MatDb')`. The original values are listed below

\begin{table}[H]
	\caption{Material properties {m}`d_rail('nmap.MatDb')`.}
	\label{tab:schemes}
	\centering
	\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}ccccc}
	    \hline
		 		&Wheel/rail	& Bracket/bolt	&Pad	& Sleeper \\
		 \hline
	 Material  &Steel &Steel &Rubber &Concrete  \\    
      Notation & UIC30 & E120/B41 & Pad& S10 (sleeper)\\    
      Module de Young ($GPa$) & 200 &200&0.014 &20\\    
      Coefficient de poisson)& 0.3 &0.3&0.45& 0.2\\   
      Density ($kg/m^3$) & 7829&7829&1000 &2000\\
      Behavior & Elastic&Elastic&Viscoelastic &Elastic\\
		\hline
	\end{tabular*}
\end{table}

Entre deux tronçons de rail, au niveau de l'éclisse, il y a un espace vide que l'on appelle \textbf{le gap}. Pour combler ce trou, nous utilisons des joints qui se différencient par leurs propriétés physiques mais aussi par leur épaisseur. En annexe, nous trouvons un tableau avec toutes les données nécessaires sur les différents types de joints.

xxx lower position of nodes within gap xxx 

\cssection{Boundary conditions}{bc}

Bondary conditions are 
\begin{enumerate}
\item {m}`sym} : half track symmetry {s}`y==717.5 -DOF 2` as y=0 corresponds to xxx middle of rail foot ? 
\end{enumerate}

\cssection{Sample meshing details}{sample}

La première étape consiste à créer des sections des différentes pièces présentées sur la figure \ref{sections} avant d'être extrudées à l’aide d’une longueur caractéristique. Elles ont été construites à partir des dimensions exactes pour se rapprocher davantage de la réalité. Pour que les nœuds coïncident parfaitement entre les pièces comme par exemple le rail et les éclisses, des petites zones de transitions ont été mises en place de manière visible sur les extrémités des sections. \\ L’assemblage de l’éclisse avec le rail se fait à l’aide de boulons permettant de fixer ces deux pièces entre elles. La modélisation de ces boulons est primordiale et ne peut pas être négligée. Ces derniers sont des zones volumiques contrairement aux autres et ainsi pour obtenir une pièce initiale incluant les boulons, nous devons faire appel aux autres sections en les extrudant afin créer une pièce volumique illustrée sur la figure \ref{sections}. De plus, il est nécessaire de créer un trou dans l’âme du rail et des éclisses afin de permettre l'assemblage des boulons.\\\\ L'étape suivante est la construction du maillage global. Pour cela, une méthode itérative a été utilisé  extrudant les sections à l’aide de la longueur caractéristique. Les pièces sont appelées en fonction des paramètres d’entrées et assemblées les unes aux autres comme illustré sur la figure \ref{maillagefinale}. Dans la première partie de l’assemblage, nous nous intéressons pas encore à la zone de contact. 

\begin{figure}[H]
    	\centering
       \includegraphics[width=0.6\textwidth]{modelschema1}
        \caption{Les étapes de l'assemblage du maillage final}
        \label{maillagefinale}
\end{figure}

La construction est réalisée étape par étape en s’aidant des sections créées précédemment. Cette stratégie permet d'avoir une flexibilité sur la modification de la taille du maillage en jouant sur la longueur caractéristique, les dimensions de l'éclisse, l'état du rail et le nombre de boulons. L’objectif est de paramétriser le plus possible afin d’avoir une large gamme d'étude par la suite.


\cssection{Rail refinement}{refinerail}
%Mettre la raideur de contact, taille de la surface, de la zone de contact 1.5-2cm de longueur et 1 cm largeur.
% 

Une fois que la géométrie et le maillage du système ont été réalisés, il faut mettre en place la zone de contact entre la roue et le rail. Toute la complexité réside dans cette zone précise et demande un maillage beaucoup plus raffiné. L’objectif est de prendre une surface bien délimitée sur le dessus du rail et de modifier la taille du maillage en conservant la même méthodologie. Le raffinement dans les zones de contact doit être le plus limité possible pour réduire au maximum la taille du maillage ainsi que limiter la zone de recherche du contact.
%Le plan d'action est de 
%le but est de 
%l'idée est de
\begin{figure}[H]
    	\centering
       \includegraphics[width=0.7\textwidth]{modelschema2}
       
        \caption{Le maillage de la zone de contact.}
        \label{contactm}
\end{figure}
\newpage
Le raffinement de cette zone est réalisé en plusieurs étapes : 
\begin{itemize}
\item on extrait d’un maillage initial une zone bien délimitée correspondant à la zone de contact;
\item les éléments de la zone à raffiner en contact avec le maillage extérieur sont transformés en éléments de transition (cf. fig \ref{contactm});
\item les éléments strictement à l’intérieur de la zone à raffiner sont divisés en plusieurs éléments plus petits. Nous avons divisé les lignes d’un cube par 3 afin d’obtenir 9 petits cubes. 

\item profile adjustement done at the end xxx 
\end{itemize}

Nous pouvons remarquer que la taille de maillage de la zone de contact est divisée par 3 par rapport au cube dans la zone lointaine. Cependant, dans l’ancien modèle sur EUROPLEXUS, ce raffinement dans cette zone était divisé par 9. La taille de maillage de l’éclisse de l’ancien maillage correspond à la taille de maillage de la zone de contact du nouveau.  Cette nouvelle approche se justifie par l’amélioration de la mise en place du contact entre la roue et le rail.\\
Plus généralement, l’ancien modèle possédait un nombre de nœuds beaucoup plus élevé que le nouveau provenant du fait qu’ils ont opté pour un maillage raffiné à partir de la zone de début et fin de l'éclisse. Ils ont essayé de réduire le raffinement sur une zone plus restreinte et de voir si la convergence des résultats est correcte. En l'occurrence, notre étude est focalisée sur les déformations au niveau de la zone de contact. C’est pourquoi un tel maillage n’est pas optimal dans notre cas. Un raffinement fin de l’éclisse, du boulon et de la surface de contact augmente de manière conséquente le temps de calcul.  Un des objectifs a été de réduire ce temps de calcul en se focalisant seulement sur les zones les plus sensibles et visées.


Après rafinnement, la géométrie réalisée pour le profil de rail ne correspond pas exactement à la réalité.  A cause du maillage, le raffinement de départ est légérement facettisé.  Les points de contact sont sur les facettes et comme illustré sur la figure \ref{profilrail2} la tête du rail est anguleuse et donc contradictoire avec la réalité. \\ Afin de remédier à ce problème, une stratégie a été mise en place consistant à projeter un profil de rail réel sur le maillage. Ce profil est représenté par les points bleus sur la figure \ref{profilrail2}.

\begin{figure}[H]
		\centering
    	\includegraphics[width=.3\textwidth]{ligne1}
       \includegraphics[width=.3\textwidth]{ligne2}
    \caption{Le profil d'un rail sur le maillage.  \label{profilrail2}}    
\end{figure}

Dans les zones de départ le maillage grossier est proche du profil avec un écart maximum de 1/10 de millimètre. Une projection orthogonale des nœuds de la surface du maillage raffiné a été faite sur le profil réel. La conséquence de cette approximation est que la hauteur de la surface de contact est faussée  de quelques millimètres.


\cssection{Details of bracket configurations}{diffeclisse}

Au fil du temps la voie ferrée a subit des transformations et des améliorations sur l'ensemble du territoire. Avec plus de 51 217 kilomètres de voies ferrées en France, il est difficilement possible d'avoir un modèle unique d'éclisse. Elles se différencient généralement de leur taille et du nombre de boulons présents. Généralement, les éclisses possèdent 4 ou 6 boulons comme illustré sur la figure \ref{diffeclisse}. Par contre, la mise en place des traverses autour de l'éclisse peut se différencier d'une éclisse à une autre. Afin d'avoir un large choix d'études, une automatisation a été mise en place dans le code permettant ainsi de passer d'un maillage à un autre avec seulement un appel. Sur la figure \ref{diffeclisse}, nous pouvons voir les différents cas utilisés et leur modèle numérique. Par exemple, le modèle $B41$ correspond à une éclisse possédant 4 boulons avec en dessous la présence de deux traverses proches entre elles.
\begin{figure}[H]
    	\centering
       \includegraphics[width=0.5\textwidth]{diffeclisse1}
        \caption{Predefined bracket configurations, {m}`Ec41} (previously B41, E120), {m}`Ec42}, {m}`Ec61}, {m}`Ec62}}
        \label{diffeclisse}
\end{figure}
Des études pourront être menées pour voir l'influencer de la présence des traverses sous l'éclisse et d'en déduire des conclusions pertinentes.

\begin{itemize}
\item les éclisses sont en contact avec le rail avec seulement deux zones de contact et non sur toute une surface comme dans l’ancien modèle présenté dans le rapport \cite{PerniaSanchez_00}%{biblio12}; 
\item l'intercalaire collant n’est pas modélisé car les contraintes de serrage à la fin sont identiques et ce qui importe ce sont les contraintes de serrage incluses dans les contraintes des rails. 
\end{itemize}

\cssection{Pad models}{PadMesh}

{m}`PadFuSN\{io4,f0-18\}} xxx 


\cssection{Wheel models}{WheelMesh}


Implemented wheel types are selected by an URN of the type {m}`W2\{xaza,in104\}}

\begin{itemize}
\item {W0} no wheel, used for moving load or frequeny response computations 
\item {W1} wheel with node to surface contact. 
\item {W2} wheel with gauss to surface contact. 
\end{itemize}


\cssection{Contact models}{Ctc}

Master gauss points of contact elements are positionned on the wheel side with {m}`ProId=3}, to avoid using very fine meshes in explicit dynamics, the surface is refined using a fine integration rule xxx. The surface for potential matching is defined as a set of display elements with  {m}`ProId=2}

For volume/volume configurations one considers {m}`ScLd} contact (surface contact large displacement) implemented in {m}`nl\_contact}.  
%
\begin{verbatim}
struct('Kc',[-23e8],'delta',[20], 'CtcForm',[13], ...
'Fl',[1000000], ...
'xgap',[-7.5 7.5])
il=[3 fe_mat('p_contact','MM',2) -2 104 1 1 0 0 0 0 0]
\end{verbatim}


In Dynavoie, the contact is implemented in {m}`nl\_dynavoie} with a default $kc=1e9 N/m$ and a {m}`.traj.CDofPos} field containing the precomputation of DOF positions. This should be moved to the ScLd format. Verify the latest Oscar case for line search. 

For rail19, and a surface of $129 mm^2$ For Kc too small the surface increases see \ser{HertzDprs}, xxx 
