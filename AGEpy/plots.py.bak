import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec

def CellPlot(df, output_file=None, gene_expression="log2FC", figure_title="CellPlot", pvalCol="elimFisher", lowerLimit=None, upperLimit=None, colorBarType='Spectral'):
    """
    Python implementation of the CellPlot from the CellPlot package for R.
    -inf or inf enrichments will come out as min found float or max found float, respectively.

    :param df: pandas dataframe with the following columns - 'Enrichment', 'Term', and 'log2fc'.
               For log2fc each cell must contain a comma separated string with the log2fc for the genes enriched in the respective term.
               eg. '-inf,-1,2,3.4,3.66,inf'
    :param output_file: prefix for an output file. If given it will create output_file.CellPlot.svg and output_file.CellPlot.png
    :param gene_expression: label for the color gradiant bar.
    :param figure_title: Figure title.
    :param pvalCol: name of the column containing the p values to determine if the terms should be marked as NS - not significant, use None for no marking
    :param lowerLimit: lower limit for the heatmap bar (default is the 0.1 percentile)
    :param upperLimit: upper limit for the heatmap bar (default is the 0.9 percentile)
    :param colorBarType: type of heatmap, 'Spectral' is dafault, alternative eg. 'seismic'
    :returns: a matplotlib figure
    """
    limits=pd.DataFrame(df['log2fc'].str.split(",").tolist())
    limits=limits.as_matrix().flatten()
    limits=pd.DataFrame(limits)
    limits[0]=limits[0].astype(str)
    limits=limits[limits[0]!="None"][limits[0]!=""]

    try:
        limits=[float(x) for x in limits[0].tolist()]

    except ValueError,e:
        print "error",e,"on line"

    if upperLimit:
        maxFC=upperLimit
    else:
        maxFC=np.percentile(limits,90)

    if lowerLimit:
        minFC=lowerLimit
    else:
        minFC=np.percentile(limits,10)

    #maxFC=np.percentile(limits,90)
    #minFC=np.percentile(limits,10)

    cmap = matplotlib.cm.get_cmap(colorBarType)
    norm = matplotlib.colors.Normalize(vmin=minFC, vmax=maxFC)

    if len(df)>10:
        siz=len(df)*3/10
    elif len(df)==1:
        siz=1
    elif len(df)==2:
        siz=2
    else:
        siz=3


    fig = plt.figure(figsize=(8, siz))
    #fig.suptitle(figure_title, fontsize=24, fontweight='bold')

    ax1 = fig.add_axes([0.05, 3.5/( float(siz)*float(10)/float(3) ), 0.9, 2])
    ax2 = fig.add_axes([0.05, 1.5/( float(siz)*float(10)/float(3) ), 0.9, 1.5/( float(siz)*float(10)/float(3) )])
    arrangment=np.arange(len(df))+.5
    enr=df['Enrichment'].tolist()
    enr=[x for x in enr if str(x) != str(float("inf"))]
    enr=[x for x in enr if str(x) != str(float("-inf"))]

    m=max(enr)

    maxE=max(enr)
    minE=min(enr)

    def getINFs(x,maxFC=maxFC,minFC=minFC):
        if x == str(float("inf")):
            return maxFC
        elif x == str(float("-inf")):
            return minFC
        else:
            return x

    def fix_enrichment(x,minE=minE,maxE=maxE):
        if str(x) == str(float("inf")):
            return maxE
        elif str(x) == str(float("-inf")):
            return minE
        else:
            return x

    df['Enrichment']=df['Enrichment'].apply(lambda x: fix_enrichment(x))


    ax1.barh(arrangment, df['Enrichment'].tolist(), color='white', edgecolor='black')#range(0,len(test))
    for i,pos in zip(df.index.tolist(),arrangment):
        fcs=df.ix[i,'log2fc'].split(",")
        fcs=pd.DataFrame(fcs)
        fcs[0]=fcs[0].astype(str)
        fcs[0]=fcs[0].apply(lambda x: getINFs(x))

        #fcs=fcs[fcs[0]!=""].astype(float)[0].tolist()
        fcs=fcs.astype(float)[0].tolist()

        try:
            w=float(df.ix[i,'Enrichment'])/float(len(fcs))
        except:
            print df.ix[i,]
        p=0
        fcs.sort(key=float)
        for f in fcs:
            #if float(f) > maxFC:
            #    f=maxFC
            #if float(f) < minFC:
            #    f=minFC
            #ax1.barh(pos, w, left=p, color=cmap(norm(float(f))), edgecolor='black')
            ax1.barh(pos, w, left=p, color=cmap(norm(float(f))), edgecolor=cmap(norm(float(f))))
            p=p+w
        if pvalCol:
            if df.ix[i,pvalCol] < 0.05:
                barAn=len(fcs)
            else:
                barAn=str(len(fcs))+" (NS)"
        else:
            barAn=len(fcs)
        ax1.text(df.ix[i,'Enrichment']+m*.02, pos+0.25, barAn, ha='left', va='bottom')

    ax1.set_yticks(arrangment+0.4)
    ax1.set_yticklabels(df['Term'].tolist())

    ax1.tick_params(
        axis='y',
        which='both',
        left='off',
        right='off',
        labelleft='on')

    ax1.tick_params(
        axis='x',
        which='both',
        bottom='off',
        top='on',
        labelbottom='off',
        labeltop='on')

    ax1.set_ylim(ymax = max(arrangment) + 1.5 ) #1.5
    ax1.set_xlabel("GO Term Enrichment")
    ax1.xaxis.set_label_position('top')

    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)

    cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap,norm=norm, orientation='horizontal')
    if not lowerLimit:
        if not upperLimit:
            cb1.set_label(gene_expression+'\n(0.1-0.9 percentiles)\n\n\n'+figure_title)
        else:
            cb1.set_label(gene_expression+'\n\n\n'+figure_title)
    else:
        cb1.set_label(gene_expression+'\n\n\n'+figure_title)

    #plt.subplots_adjust(top=8.5)

    #ax1.annotate('a fractional title', xy=(.025, .975), xycoords='figure fraction', horizontalalignment='center', verticalalignment='top', fontsize=20)
    #plt.subplots_adjust(top=500)
    #ax1.text(2, 2, 'right bottom',
    #    horizontalalignment='center',
    #    verticalalignment='center',
    #    transform=ax1.transAxes)

    if output_file:
        plt.savefig(output_file+".CellPlot.png",dpi=300,bbox_inches='tight', pad_inches=0.1,format='png')
        plt.savefig(output_file+".CellPlot.svg",dpi=300,bbox_inches='tight', pad_inches=0.1,format='svg')

    return fig

def SymPlot(df,output_file=None,figure_title="SymPlot",pvalCol="elimFisher"):
    """
    Python implementation of the SymPlot from the CellPlot package for R.
    -inf or inf enrichments will come out as min found float or max found float, respectively.

    :param df: pandas dataframe with the following columns - 'Enrichment', 'Significant', 'Annotated', 'Term', and 'log2fc'.
               For log2fc each cell must contain a comma separated string with the log2fc for the genes enriched in the respective term.
               eg. '-inf,-1,2,3.4,3.66,inf'
    :param output_file: prefix for an output file. If given it witll create output_file.SymPlot.svg and output_file.SymPlot.png
    :param figure_title: Figure title.
    :param pvalCol: name of the column containing the p values to determine if the terms should be marked as NS - not significant, use None for no marking
    :returns: a matplotlib figure
    """
    maxAn=df['Annotated'].max()

    arrangment=np.arange(len(df))+.5

    def getINFs(x):
        if x == str(float("inf")):
            return 1
        elif x == str(float("-inf")):
            return -1
        else:
            return x

    enr=df['Enrichment'].tolist()
    enr=[x for x in enr if str(x) != str(float("inf"))]
    enr=[x for x in enr if str(x) != str(float("-inf"))]

    maxE=max(enr)
    minE=min(enr)

    def fix_enrichment(x,minE=minE,maxE=maxE):
        if str(x) == str(float("inf")):
            return maxE
        elif str(x) == str(float("-inf")):
            return minE
        else:
            return x

    df['Enrichment']=df['Enrichment'].apply(lambda x: fix_enrichment(x))

    limits=df['Enrichment'].tolist()
    maxFC=np.percentile(limits,90)
    minFC=np.percentile(limits,10)

    cmap = matplotlib.cm.get_cmap('Spectral')
    norm = matplotlib.colors.Normalize(vmin=minFC, vmax=maxFC)

    if len(df) >= 5:
        siz=len(df)*4/10
    else:
        size=5*4/10


    fig = plt.figure(figsize=(8, siz))
    #fig.suptitle(figure_title, fontsize=24, fontweight='bold')
    gs = gridspec.GridSpec(1, 3, width_ratios=[2,0.75,2])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    ax4 = fig.add_axes([0.125, 0.11/100*len(df), 0.775, 0.075*10/len(df)])#/100.0075

    alldown=[]
    allup=[]

    for i,pos in zip(df.index.tolist(),arrangment):
        f=df.ix[i,'Enrichment']
        sigN=df.ix[i,'Significant']
        ann=float(df.ix[i,'Annotated'])

        if ann!=maxAn:
            p=float(maxAn-ann)/2
        else:
            p=0
        ax2.barh(pos, ann, left=p, color=cmap(norm(float(f))),edgecolor=cmap(norm(float(f))))#

        fcs=df.ix[i,'log2fc'].split(",")
        fcs=pd.DataFrame(fcs)
        fcs[0]=fcs[0].astype(str)
        fcs[0]=fcs[0].apply(lambda x: getINFs(x))
        #fcs=fcs[fcs[0]!=""].astype(float)
        fcs=fcs.astype(float)
        down=len(fcs[fcs[0]<0])/ann*100
        up=len(fcs[fcs[0]>0])/ann*100
        alldown.append(down)
        allup.append(up)

        ax1.barh(pos, down, color="blue",edgecolor='blue')
        ax3.barh(pos, up, color="red",edgecolor='red')

    ax1.spines['top'].set_visible(True)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    ax1.tick_params(axis='x',which='both',bottom='off', top='on',labelbottom='off',labeltop='on')
    ax1.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')

    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    ax2.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off',labeltop='off')
    ax2.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')

    ax3.spines['top'].set_visible(True)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    ax3.tick_params(axis='x',which='both',bottom='off',top='on',labelbottom='off',labeltop='on')
    ax3.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')

    fa=10*0.1/len(df)+1
    fb=10*0.08/len(df)+1

    ax1.set_title('Downregulated (%)',y=fa)#
    ax2.set_title('Annotated\n(max=%s)' %str(maxAn),y=fb)#
    ax3.set_title('Upregulated (%)',y=fa)

    ax1.set_xlim(max(max(alldown),max(allup)), 0)
    ax2.set_xlim(0, maxAn)
    ax3.set_xlim(0, max(max(alldown),max(allup)))

    ax1.set_ylim(ymax = max(arrangment)+1.5)
    ax2.set_ylim(ymax = max(arrangment)+1.5)
    ax3.set_ylim(ymax = max(arrangment)+1.5)


    ax1.set_yticks(arrangment+0.4)
    def get_label_with_sig (df):
        termLabel=df['Term']
        if pvalCol:
            pvalue=df[pvalCol]
            if pvalue > 0.05:
                return "(NS) "+termLabel
            else:
                return termLabel
        else:
            return termLabel

    df['newLabels']=df.apply(get_label_with_sig, axis=1)

    ax1.set_yticklabels(df['newLabels'].tolist())

    cb1 = matplotlib.colorbar.ColorbarBase(ax4, cmap=cmap,norm=norm, orientation='horizontal')
    cb1.set_label('GO Term Enrichment (0.1-0.9 percentiles)\n\n\n'+figure_title)

    fig.subplots_adjust(wspace=0)

    if output_file:
        plt.savefig(output_file+".SymPlot.png",dpi=300,bbox_inches='tight', pad_inches=0.1,format='png')
        plt.savefig(output_file+".SymPlot.svg",dpi=300,bbox_inches='tight', pad_inches=0.1,format='svg')

    return fig

def NormInt(df,sampleA,sampleB):
    """
    Normalizes intensities of a gene in two samples

    :param df: dataframe output of GetData()
    :param sampleA: column header of sample A
    :param sampleB: column header of sample B

    :returns: normalized intensities
    """

    c1=df[sampleA]
    c2=df[sampleB]
    return np.log10(np.sqrt(c1*c2))

def MA(df,title,figName,c, daType="counts",nbins=10,perc=.5,deg=3,eq=True,splines=True,spec=None,Targets=None,ylim=None,sizeRed=8):
    """
    Plots an MA like plot GetData() outputs.

    :param df: dataframe output of GetData()
    :param title: plot title, 'Genes' or 'Transcripts'
    :param figName: /path/to/saved/figure/prefix
    :param c: pair of samples to be plotted in list format
    :param daType: data type, ie. 'counts' or 'FPKM'
    :param nbins: number of bins on normalized intensities to fit the splines
    :param per: log2(fold change) percentil to which the splines will be fitted
    :param deg: degress of freedom used to fit the splines
    :param eq: if true assumes for each bin that the lower and upper values are equally distant to 0, taking the smaller distance for both
    :param spec: list of ids to be highlighted
    :param Targets: list of ids that will be highlighted if outside of the fitted splines
    :param ylim: a list of limits to apply on the y-axis of the plot
    :param sizeRed: size of the highlight marker

    :returns df_: a Pandas dataframe similar to the GetData() output with normalized intensities and spline outbounds rows marked as 1.
    :returns red: list of ids that are highlighted
    """

    df_=df[df[c[0]]>0]
    df_=df_[df_[c[1]]>0]

    df_["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) ]=df_.apply(NormInt, args=(c[0],c[1],), axis=1)

    if daType=="counts":
        lowLim=np.log10(np.sqrt(10))
    elif daType=="FPKM":
        lowLim=np.log10(0.1)

    df_b=df_[df_["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) ]>lowLim ]
    df_b.reset_index(inplace=True, drop=True)

    Xdata=df_["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) ].tolist()
    Ydata=df_["log2(%s/%s)" %( str(c[1]), str(c[0]) )].tolist()

    minX=min(Xdata)
    maxX=max(Xdata)

    minX_=min(df_b["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) ].tolist())
    maxX_=max(df_b["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) ].tolist())

    df_b["bin"]=pd.cut(df_b["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) ].tolist(), nbins,labels=False)

    spl=[]
    for b in set( df_b["bin"].tolist() ):
        tmp=df_b[df_b["bin"]==b]
        Xbin = tmp["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) ].tolist()
        Xval = np.mean([max(Xbin),min(Xbin)])
        Ybin = tmp["log2(%s/%s)" %( str(c[1]), str(c[0]) )].tolist()
        YvalP=np.percentile(Ybin,100.00-float(perc))
        YvalM=np.percentile(Ybin,float(perc))
        spl.append([Xval,YvalP,YvalM])

    spl=pd.DataFrame( spl,columns=["X","Upper","Lower"],index=range(len(spl)) )

    def CheckMin(df):
        U=abs(df["Upper"])
        L=abs(df["Lower"])
        return min([U,L])

    spl["min"]=spl.apply(CheckMin, axis=1)

    coeffsUpper = np.polyfit(spl["X"].tolist(), spl["Upper"].tolist(), deg)
    coeffsLower = np.polyfit(spl["X"].tolist(), spl["Lower"].tolist(), deg)

    Xspl = np.array(np.linspace(minX, maxX, 10*nbins))

    if eq:
        coeffsUpper = np.polyfit(spl["X"].tolist(), spl["min"].tolist(), deg)
        coeffsLower = np.polyfit(spl["X"].tolist(), [ ss*-1 for ss in spl["min"].tolist()] , deg)
        YsplUpper = np.polyval(coeffsUpper, Xspl)
        YsplLower = np.polyval(coeffsLower, Xspl)

    else:
        coeffsUpper = np.polyfit(spl["X"].tolist(), spl["Upper"].tolist(), deg)
        coeffsLower = np.polyfit(spl["X"].tolist(), spl["Lower"].tolist(), deg)
        YsplUpper = np.polyval(coeffsUpper, Xspl)
        YsplLower = np.polyval(coeffsLower, Xspl)

    def checkOutbounds(df,Xspl=Xspl,coeffsUpper=coeffsUpper,coeffsLower=coeffsLower,c=c):
        x=df["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) )]
        y=df["log2(%s/%s)" %( str(c[1]), str(c[0]) )]
        if y < 0:
            v=np.polyval(coeffsLower, x)
            if y < v:
                return 1
            else:
                return 0
        else:
            v=np.polyval(coeffsUpper, x)
            if y > v:
                return 1
            else:
                return 0

    df_["OutBounds"]=df_.apply(checkOutbounds,axis=1)

    if Targets:
        if title == "Transcripts":
            red=df_[df_["OutBounds"]==1][df_["transcript_id"].isin(Targets)]["transcript_id"].tolist()
            Xdata_=df_[df_["OutBounds"]==1][df_["transcript_id"].isin(Targets)]["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) ].tolist()
            Ydata_=df_[df_["OutBounds"]==1][df_["transcript_id"].isin(Targets)]["log2(%s/%s)" %( str(c[1]), str(c[0]) ) ].tolist()
        elif title == "Genes":
            red=df_[df_["OutBounds"]==1][df_["gene_id"].isin(Targets)]["gene_id"].tolist()
            Xdata_=df_[df_["OutBounds"]==1][df_["gene_id"].isin(Targets)]["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) )].tolist()
            Ydata_=df_[df_["OutBounds"]==1][df_["gene_id"].isin(Targets)]["log2(%s/%s)" %( str(c[1]), str(c[0]) )].tolist()
    elif spec:
        if title == "Transcripts":
            red=df_[df_["transcript_id"].isin(spec)]["transcript_id"].tolist()
            Xdata_=df_[df_["transcript_id"].isin(spec)]["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) ].tolist()
            Ydata_=df_[df_["transcript_id"].isin(spec)]["log2(%s/%s)" %( str(c[1]), str(c[0]) ) ].tolist()
        elif title == "Genes":
            red=df_[df_["gene_id"].isin(spec)]["gene_id"].tolist()
            Xdata_=df_[df_["gene_id"].isin(spec)]["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) )].tolist()
            Ydata_=df_[df_["gene_id"].isin(spec)]["log2(%s/%s)" %( str(c[1]), str(c[0]) )].tolist()
    else:
        Xdata_=df_[df_["OutBounds"]==1]["normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) ].tolist()
        Ydata_=df_[df_["OutBounds"]==1]["log2(%s/%s)" %( str(c[1]), str(c[0]) ) ].tolist()
        if title == "Transcripts":
            red=df_[df_["OutBounds"]==1]["transcript_id"].tolist()
        elif title == "Genes":
            red=df_[df_["OutBounds"]==1]["gene_id"].tolist()

    fig = plt.gcf()
    fig.set_size_inches(6, 6)
    plt.scatter(Xdata,Ydata, s=2)
    plt.scatter(Xdata_,Ydata_,s=sizeRed, c='r')
    if splines:
        plt.plot(Xspl,YsplUpper, "-",lw=0.5, c='g')
        plt.plot(Xspl,YsplLower,"-", lw=0.5,c='g')

    plt.xlabel("normalized intensities (%s vs. %s)" %( str(c[0]), str(c[1]) ) )
    plt.ylabel("log2(%s/%s)" %( str(c[1]), str(c[0]) ))

    if ylim:
        plt.ylim(ylim[0],ylim[1])
    else:
        ylims=max([abs(min(Ydata)), abs(max(Ydata)) ])
        plt.ylim(-ylims*1.1,ylims*1.1)

    plt.title(title)
    plt.savefig(figName+".png",dpi=300,bbox_inches='tight', pad_inches=0.1,format='png')
    plt.savefig(figName+".svg",dpi=300,bbox_inches='tight', pad_inches=0.1,format='svg')
    plt.show()

    return df_,red
