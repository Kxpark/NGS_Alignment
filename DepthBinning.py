import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--input', '-i', required=True, help='input SAMtools depth file')
parser.add_argument('--bin', '-b', required=True, help='The bin size in bp (integer)',type=int)
parser.add_argument('--prefix','-p',required=True,help='Put a prefix, including a path, of the generated output files',type=str)
parser.add_argument('--Ymax',required=False,help='The Y max level on the output of the graph',type=int,default=30)
args = parser.parse_args()


def getBinTotal(file,bin,output,ymax):
    binnData=[]
    df= pd.read_csv(file,delimiter='\t')
    U=range(0,(max(df.POS)),bin)
    print(max(df.POS))
    numbins=int(max(df.POS)/bin)
    if len(df.columns) == 3:
        print('Only One Line')
        n=0
        for bn in range(0,(numbins-1),1):
                lowbin=U[0+n]
                upbin=U[1+n]
                binnData.append(sum(df[df.columns[2]][df['POS'].between(lowbin,upbin)])/bin)
                n=n+1
        X=list(U)
        del X[0]
        del X[(numbins-1)]
        Y=list(binnData)
        df2=pd.DataFrame({'X':X,'Y':Y})
        df2.to_csv(output+'Binned.txt',sep='\t')
        ax=df2.plot(kind='scatter',x='X',y='Y',color='b')
        ax.set_ylim(0,ymax)
        ax.set_xlabel('Base Pairs in bins of '+str(bin)+'bp')
        ax.set_ylabel('Coverage')
        ax.set_title('Depth of Reads per Bin Region')
        fig=ax.get_figure()
        fig.savefig(output+"DepthPlot.png",dpi=300)
        
    else:
        totalData={}
        lines=(len(df.columns)-2)
        print('Total Lines'+str(lines))
        for lin in range(0,lines):
             binnData=[]
             n=0
             for bn in range(0,(numbins-1),1):
                lowbin=U[0+n]
                upbin=U[1+n]
                binnData.append(sum(df[df.columns[(lin+2)]][df['POS'].between(lowbin,upbin)])/bin)
                n=n+1
             totalData[df.columns[lin+2]]=binnData
        
        keys=list(totalData.keys())
        X=list(U)
        del X[0]
        del X[(numbins-1)]
        totalData['X']=X
        df2=pd.DataFrame(totalData)
        df2.to_csv(output+'Binned.txt',sep='\t')
        ax=df2.plot(kind='scatter',x='X',y=keys[0],color='b',label=keys[0])
        for line in range(1,lines):
             df2.plot(kind='scatter',x='X',y=keys[line],ax=ax,color='r',label=keys[line])
        ax.set_ylim(0,ymax)
        ax.set_xlabel('Base Pairs in bins of '+str(bin)+'bp')
        ax.set_ylabel('Coverage')
        ax.set_title('Depth of Reads per Bin Region')
        fig=ax.get_figure()
        fig.savefig(output+"DepthPlot.png",dpi=300)
        df2.to_csv(output+'Binned.txt',sep='\t')



getBinTotal(args.input,args.bin,args.prefix,args.Ymax)
print(f'wrote {args.input} to binned file using {args.bin} bps per bin! and saved to {args.prefix}.Binnedfile.txt') 
        
