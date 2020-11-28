**Workflow for creating probtrack os2t output networks in Gephi**

Step 1)  **Creating Input Nodes and Edges Files for Gephi**
1.  First step is to create the nodes table, All nodes tables require at
    minimum 2 columns: “Id” and “Label”.
    
-   **Id:** Numbered list from 1number of regions, this includes both
    the targets and the seeds

-   **Label:** Label corresponding to that Id (i.e. row 1: Id=1,
    Label=V1 … etc). If looking specifically at the connectivity of
    seeds—target, place the seeds at the bottom of the list (i.e. 180
    regions in glasser atlas targets, then seed labels begin at
    Id=181, Label=Seed\#1Name)

-   You may add additional columns in the Nodes table that describe some
    property of that Id,label pair. (i.e. common to add “Modularity
    Class” that describes clusters that targets belong two, with
    additional classes for each seed (i.e. Id=2,3,4 all have Modularity
    Class 2 … you may also add a column of labels for these Modularity
    Classes as Modularity Class Label (i.e. Id=2,3,4 all have Modularity
    Class = 2 and Modularity Class Label = Early Visual Cortex), along
    with the accompanying “**color”** column that denotes RBG values for
    this modularity class. Any other descriptors may be added to the
    nodes table

-   Save the Nodes table as a **csv**

-   Example Nodes Table:
```
  Id   Label   color         Modularity Class
  ---- ------- ------------- ------------------
  1    V1      180,182,255   1
  2    V2      0,115,155     2
  3    V3      0,115,155     2
  4    V4      0,115,155     2
  5    V6      173,192,0     3
  6    V3B     173,192,0     3
  7    V6A     173,192,0     3
  8    IPS1    173,192,0     3
  9    V7      173,192,0     3
  10   V3A     173,192,0     3
```

2.  Creating the edges table, All edges tables require at minimum 2
    columns: “Source” and “Target”. Note that since there are 3
    different seeds (Cluster01, Cluster02, Cluster03) that each have 180
    targets (Glasser atlas targets 1-180), each column length is
    180\*number of seeds (i.e. 180x3) with the first 180 entries
    corresponding to Cluster01, etc.

-   **Source:** Id corresponding to the source node (i.e. 181)

    -   Best practice is to also include the next column as **SLabel**:
        Label corresponding to that source node (I.e. Id 181 = Cluster01

-   **Target:** Id corresponding to the target node (i.e. 1)

    -   Best practice is to also include the next column as **TLabel:**
        Label corresponding to that target node (i.e. Id 1 = V1)

-   You may add any additional columns in the Edges table that describe
    some property of the Source—Target pair (i.e. ConnectivityScore,
    ConnectivityScoreNormalized, RankofConnectivityScore). You can use
    these later to filter your graphs in Gephi.

-   A non-essential but important additional column to include is
    **Weight**: this column defines the edge weight for that
    Source-Target pair. Note that weight must be &gt; 0, as weight = 0
    will be considered to be a non-existent edge. In diffusion, this
    weight column could be the for example the normalized connectivity
    score scaled from 1100. In my case, Weight =
    ConnectivityScoreNormalized (per seed) \* Rank\_in\_Cluster (per
    seed), rounded to nearest integer, then adding 1 to ensure all
    weights are &gt; 0

-   Save the Edges table as a **csv**

-   Example Edges table (note that I sorted this table to view the
    highest weight at the top, just for simplicity in viewing this
    example):
```
  Source   SLabel      Target   TLabel   ConnectivityScore   ConnectivityScoreNormalized   Rank\_in\_Cluster   Weight
  -------- ----------- -------- -------- ------------------- ----------------------------- ------------------- --------
  181      Cluster01   70       AVI      4400.303            0.486185                      180                 89
  182      Cluster02   70       AVI      1737.989            0.222274                      180                 41
  183      Cluster03   64       STGa     337.4957            0.183807                      180                 34
  182      Cluster02   73       PoI1     1400.25             0.17908                       179                 33
  183      Cluster03   73       PoI1     300.3282            0.163565                      179                 30
  181      Cluster01   73       PoI1     1358.506            0.1501                        179                 28
  183      Cluster03   70       AVI      268.062             0.145992                      178                 27
  183      Cluster03   80       H        249.5958            0.135935                      177                 25
  182      Cluster02   89       TE1a     1012.787            0.129527                      178                 24
  181      Cluster01   89       TE1a     950.5206            0.105022                      178                 20
```

Step 2)  **Loading the Nodes and Edges into Gephi**

1.  Within Gephi, open a new project and navigate to the “Data
    Laboratory” tab.

2.  Select “Import Spreadsheet” and select the NodesTable.csv
    you created. Ensure separator is set to comma, Import as: is set to
    Nodes table, and charset to UTF-8. Click next and ensure the columns
    you want to import are included and their datatype is selected
    accordingly (specifically for “color”, this should be of
    type String). Click next and ensure that graph type is “Mixed”,
    click more options and ensure Auto-scale and self-loops are
    selected, while Create missing nodes should be deselected, Edges
    merge strategy should be “Sum”. Finally, ensure “Append to existing
    workspace” is selected and click ok. You will now see your important
    table in the data laboratory “Nodes” tab.

3.  Select “Import Spreadsheet” and select the EdgesTable.csv
    you created. Ensure separator is set to comma, Import as: is set to
    Edges table, and charset to UTF-8. Click next and ensure the columns
    you want to import are included and their datatype is
    selected accordingly. Click next and sure that graph type is either
    “Directed or Undirected”… Directed is if the sourcetarget is
    unidirectional, and Undirected if source—target bidirectionally is
    valid, click more options and ensure Auto-scale and self-loops are
    selected, while Create missing nodes should be deselected, Edges
    merge strategy should be “Sum”. Finally, ensure “Append to existing
    workspace” is selected and click ok. You will now see your important
    table in the data laboratory “Edges” tab.

Step 3)  **Viewing and Manipulating the Network in Gephi**

1.  At this point, you can move to the “Overview” tab and visualize the
    network with its default settings. As viewed in the example below:

    ![](media/image1.png){width="5.504201662292213in"
    height="2.964388670166229in"}

2.  In the toolbar at the bottom, click the Dark “T” to show the
    node labels. The leftmost slider bar controls viewing size of the
    edge weights, while rightmost slider bar controls node label size.

3.  Although chord plots can be created with ForceAtlas by setting the
    force to something high e.g. 1000 or 10000 (depending on data size)
    and attraction to a 1/1000 \* force (e.g 1 or 10). **I recommend
    going into Tools&gt;Plugins&gt;Available Plugins, and installing
    Circular Layout and Dual Circle Layout by Matt Groeninger.**

4.  Using Dual Circle Layout (selected from Layout tab within the
    Overview tab): To create a chord graph with all nodes on the outside
    of the circle, except for the Ids at the end of nodes table
    (corresponding to 3 seeds in my case). Select Upper Order Count = 3,
    and Order Nodes by “Node ID”. Example of the result provided below:

    ![](media/image2.png){width="5.537814960629921in"
    height="2.959417104111986in"}

5.  Now that the primary graph chord has been created, it needs to be
    made more readable.

a)  The first thing we want to do is color the nodes based on the RGB
    values provided in the NodesTable.csv class. **To do this, you must
    have the Give Colors To Nodes plugin installed (click
    Tools&gt;Plugins&gt;Available Plugins and install Give Colors
    To Nodes)**. After a restart you will see a rainbow colorwheel at
    the bottom of the graph toolbar in the overview panel. Simply click
    this icon and the nodes will be colored automatically.

    ![](media/image3.png){width="5.346882108486439in"
    height="2.8739501312335958in"}

b)  The node size may also be changed by selecting the second icon on
    the right side of this tab (interweaved circles), selecting ranking,
    and selecting Degree as this attribute (the degree for each of the
    nodes on the outside is 3, while the nodes on the inside each have a
    degree of 180. Play around with different Min and Max sizes, I have
    used Min=10 and Max=20 in the below example. The label size may also
    be changed by selecting the fourth icon on the right side of this
    tab (TT), selecting ranking, and selecting Degree as this attribute.
    I have used Min=1 and Max=1.2 in the below example.

    ![](media/image4.png){width="5.3697473753280835in"
    height="2.880501968503937in"}

Step 4) **Filtering the Graph:** 
    
    Right now it is somewhat difficult to
    visually see the most important edge weights since each cluster in
    the center of the circle has 180 edges shown. If for example we only
    want to view the top 5 ranked edges for each cluster, we can filter
    the graph… In the filters tab within the overview tab, select
    attributesrange, and double click on the attribute you want to
    filter on (in this case we want to filter on rank\_in\_cluster). You
    will see this filter now in the Queries tab below. Since we only
    want the top 5 for each cluster, we set the lower bound to 176 (the
    filter is inclusive, so this will include 176,177,178,179,180) and
    then click select and filter. The filter is now continuously applied
    to this data. The two images below illustrate the parameters to set
    when filtering, and what the resulting new workspace looks like:

    ![](media/image5.png){width="5.68669072615923in"
    height="3.0420166229221346in"}

    ![](media/image6.png){width="5.706469816272966in"
    height="3.0672265966754155in"}

-   You will node that the edge colours are actually a combination of
    the nodes they connect. Don’t worry, this will be resolved when
    viewing the final graph in the preview tab.

Step 5)  **Saving the final graph**
1.  Select the preview tab while the filter is being continuously
    applied

2.  Under Nodel Labels: Select Show Labels, set font size to Arial 5
    Bold (or something similar to ensure that the node labels do not
    overlap), and make sure proportional size is selected.

3.  Under Edges: Select Rescale weight, and change the Min and Max
    rescaled weight to something readable (in this example, I have
    rescaled from Min=1, Max=20). Additionally, select the ellipses in
    “Color” and select “Source”, this changes the edge colours to
    the source. You may also deselect curved if you prefer straight
    edges

4.  Click Refresh at the bottom. If you are happy with the result:

5.  Click Export: SVG/PDF/PNG and save as a **pdf**

![](media/image7.png){width="6.5in" height="3.4972222222222222in"}

Step 6)  **Manually Adding in Modularity Class labels:**

-   Unfortunately, Gephi does not have a good way to manually add in the
    overarching modularity class labels while still maintaining the
    inner node labels, so this has to be manually done in powerpoint at
    the moment. Although once this has been done once, any additional
    graphs created can just be placed within the bounds of this
    template, as long as they were saved with the same settings
    in gephi.

-   This yields the final output below:

    ![](media/image8.png){width="5.915043744531934in"
    height="5.8632239720034995in"}


