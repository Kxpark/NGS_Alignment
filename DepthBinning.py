import pandas as pd
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt

parser = ArgumentParser()
parser.add_argument('--input', '-i', required=True, help='input SAMtools depth file')
parser.add_argument('--bin', '-b', required=True, help='The bin size in bp (integer)',type=int)
parser.add_argument('--prefix','-p',required=True,help='Put a prefix, including a path, of the generated output files',type=str)
args = parser.parse_args()


def getBinTotal(file,bin,output):
    binnData=[]
    df= pd.read_csv(file,delimiter='\t')
    U=range(0,(max(df.POS)),bin)
    print('Chromosome Length '+ str(max(df.POS)))
    numbins=int(max(df.POS)/bin)
    if len(df.columns) == 3:
            print('Only One Line')
            n=0
            for bn in range(0,(numbins-1),1):
                    lowbin=U[0+n]
                    upbin=U[1+n]
                    binnData.append(sum(df[df.columns[2]][df['POS'].between(lowbin,upbin)])/bin)
                    n=n+1
            X=np.array([U])
            X=np.delete(X,[0,numbins])
            Y=np.array(binnData)
            ymax=((np.std(Y)*5)+(np.mean(Y)))
            plt.ylim(0,ymax)
            plt.scatter(X,Y)
            plt.xlabel('Base Pairs in bins of '+str(bin)+'bp')
            plt.ylabel('Coverage')
            plt.title('Depth of Reads per Bin Region')
            plt.savefig(output+"DepthPlot.png",dpi=300)
            pd.DataFrame({'Bin':X,'Reads':Y}).to_csv(output+'Binned.txt',sep='\t')
    else:
        totalData={}
        lines=(len(df.columns)-2)
        print('Total Lines'+str(lines))
        for lines in range(0,lines):
            binnData=[]
            n=0
            for bn in range(0,(numbins-1),1):
                lowbin=U[0+n]
                upbin=U[1+n]
                binnData.append(sum(df[df.columns[2]][df['POS'].between(lowbin,upbin)])/bin)
                n=n+1
                totalData[df.columns[lines+2]]=binnData
            
        keys=list(totalData.keys())
        X=np.array([U])
        X=np.delete(X,[0,numbins])
        X=list(X)
        totalData['X']=X
        df2=pd.DataFrame(totalData)

        ax=df2.plot(kind='scatter',x='X',y=keys[0],color='b',label=keys[0])
        for lines in range(0,lines):
            df2.plot(kind='scatter',x='X',y=keys[lines],ax=ax,color='r',label=keys[lines])
        ax.set_xlabel('Base Pairs in bins of '+str(bin)+'bp')
        ax.set_ylabel('Coverage')
        ax.set_title('Depth of Reads per Bin Region')
        plt.savefig(output+"DepthPlot.png",dpi=300)
        df2.to_csv('Binned.txt',sep='\t')






getBinTotal(args.input,args.bin,args.prefix)
print(f'wrote {args.input} to binned file using {args.bin} bps per bin! and saved to {args.prefix}.Binnedfile.txt') 
        