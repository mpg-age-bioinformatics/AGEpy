## ___CellPlot___

Python implementation of the CellPlot from the CellPlot package for R.
-inf or inf enrichments will come out as min found float or max found float, respectively.

**`CellPlot(df, output_file=None, gene_expression="log2FC", figure_title="CellPlot", pvalCol="elimFisher", lowerLimit=None, upperLimit=None, colorBarType='Spectral')`**

* **`df`** pandas dataframe with the following columns - 'Enrichment', 'Term', and 'log2fc'.
           For log2fc each cell must contain a comma separated string with the log2fc for the genes enriched in the respective term.
           eg. '-inf,-1,2,3.4,3.66,inf'
* **`output_file`** prefix for an output file. If given it will create output_file.CellPlot.svg and output_file.CellPlot.png
* **`gene_expression`** label for the color gradiant bar.
* **`figure_title`** Figure title.
* **`pvalCol`** name of the column containing the p values to determine if the terms should be marked as NS - not significant, use None for no marking
* **`lowerLimit`** lower limit for the heatmap bar (default is the 0.1 percentile)
* **`upperLimit`** upper limit for the heatmap bar (default is the 0.9 percentile)
* **`colorBarType`** type of heatmap, 'Spectral' is default, alternative eg. 'seismic'
* **`returns`** a matplotlib figure

```python
>>> import AGEpy as age
>>> print df.head()

Term  Annotated  Enrichment  \
0          GO:0008544~epidermis development         38    4.006021   
1               GO:0043588~skin development         33    4.359840   
2         GO:0045087~innate immune response         61    2.385984   
3               GO:0006952~defense response         90    1.913315   
4  GO:0009605~response to external stimulus        113    1.736641   

ease                                             log2fc  
0  1.193931e-12  1.13845,0.771811,0.926561,0.578588,-0.694105,1...  
1  4.757460e-12  1.13845,0.926561,-0.694105,1.48945,0.94486,-1....  
2  5.609421e-10  -1.91507,-0.630414,-1.87466,-0.898252,0.458041...  
3  2.238959e-09  -0.538926,0.667335,-1.91507,-0.630414,-1.87466...  
4  3.051460e-09  0.667335,-1.91507,-0.630414,1.46227,0.755911,-...  

>>> cellplot=age.CellPlot(df[:20], "cellplot",  "log2(mt/wt)", "mutant 1", \
pvalCol="ease", colorBarType="bwr", lowerLimit=-1.25,upperLimit=1.25)
```
![cellpot](cellplot.CellPlot.png)
___

## ___SymPlot___

Python implementation of the SymPlot from the CellPlot package for R.
-inf or inf enrichments will come out as min found float or max found float, respectively.

**`SymPlot(df,output_file=None,figure_title="SymPlot",pvalCol="elimFisher")`**

* **`df`** pandas dataframe with the following columns - 'Enrichment', 'Significant', 'Annotated', 'Term', and 'log2fc'. 'Annotated'i stands for number of genes annotated with the respective GO term. As reported in DAVID by listHits. For log2fc each cell must contain a comma separated string with the log2fc for the genes enriched in the respective term. eg. '-inf,-1,2,3.4,3.66,inf'
* **`output_file`** prefix for an output file. If given it witll create output_file.SymPlot.svg and output_file.SymPlot.png
* **`figure_title`** Figure title.
* **`pvalCol`** name of the column containing the p values to determine if the terms should be marked as NS - not significant, use None for no marking
* **`returns`** a matplotlib figure

```python
>>> import AGEpy as age
>>> symplot=age.SymPlot(df[:20],"symplot", "mutant 1",pvalCol="ease")
```
![sympot](symplot.SymPlot.png)
___
## ___MA___

Plots an MA like plot.

**`MA(df, title, figName, c, daType="counts", nbins=10, perc=.5, deg=3, eq=True, splines=True, spec=None, Targets=None, ylim=None, sizeRed=8)`**

* **`df`** dataframe output of GetData()
* **`title`** plot title, 'Genes' or 'Transcripts'
* **`figName`** /path/to/saved/figure/prefix
* **`c`** pair of samples to be plotted in list format
* **`daType`** data type, ie. 'counts' or 'FPKM'
* **`nbins`** number of bins on normalized intensities to fit the splines
* **`per`** log2(fold change) percentil to which the splines will be fitted
* **`deg`** degress of freedom used to fit the splines
* **`eq`** if true assumes for each bin that the lower and upper values are equally distant to 0, taking the smaller distance for both
* **`splines`** plot splines, default=True
* **`spec`** list of ids to be highlighted
* **`Targets`** list of ids that will be highlighted if outside of the fitted splines
* **`ylim`** a list of limits to apply on the y-axis of the plot
* **`sizeRed`** size of the highlight marker
* **`returns df_`** a Pandas dataframe similar to the GetData() output with normalized intensities and spline outbounds rows marked as 1.
* **`returns red`** list of ids that are highlighted

```python
>>> import AGEpy as age
>>> print df.head()

gene_id                    gene  wt0  wt20  log2(wt20/wt0)  \
0  ENSG00000223972                 DDX11L1  0.0   0.0             NaN   
1  ENSG00000243485  MIR1302-2,RP11-34P13.3  0.0   0.0             NaN   
2  ENSG00000274890  MIR1302-2,RP11-34P13.3  0.0   0.0             NaN   
3  ENSG00000268020                  OR4G4P  0.0   0.0             NaN   
4  ENSG00000240361                 OR4G11P  0.0   0.0             NaN   

p_value  q_value significant  
0      1.0      1.0          no  
1      1.0      1.0          no  
2      1.0      1.0          no  
3      1.0      1.0          no  
4      1.0      1.0          no  

>>> madf1,sig1=age.MA(dge_, 'Genes',"MA1",["wt0","wt20"], daType="FPKM")
```
![ma1](MA1.png)
```python
>>> sigGenes=df[df["significant"=="yes"]]["gene_id"].tolist()
>>> madf2,sig2=age.MA(dge_, 'Genes',"MA2", ["wt0","wt20"], splines=False, daType="FPKM",spec=sigGenes)
```
![ma2](MA2.png)
```python
>>> madf3,sig3=age.MA(dge_, 'Genes',"MA3", ["wt0","wt20"], splines=True, daType="FPKM",Targets=sigGenes)
```
![ma3](MA3.png)
___
