
import numpy as np
dataset_list=[]
with open("../data/data_list","r") as whitelist:
    for line in whitelist.readlines():
        dataset=line.split(" ")[0]
        dataset_list.append(dataset)
    whitelist.close()
for para in ["","_para"]:

    perfTable = [["dataset", "ftm", "non progressive", "progressive"]]
    perfTableCP = [["dataset", "ftm", "non progressive", "progressive"]]
    for i,dataset in enumerate(dataset_list):
        dataset=dataset.split(".")[0].split("\n")[0]
        print(dataset)
        if(dataset=="ctBones"):
            perfTable.append(["foot"])
        elif(dataset=="at"):
            perfTable.append(["AT"])
        elif(dataset=="ethanediol"):
            perfTable.append(["ethaneDiol"])
        else:
            perfTable.append([dataset])
        perfTableCP.append([dataset])

        # FTM PD
        timeFTMSum=0
        write=1
        values=[]
        for itry in range(0,12):
            timeFTM=0
            count=0
            with open("./outputs_pd/ftm_"+dataset+para+"_try"+str(itry), "r") as f:
                for line in f.readlines():
                    line=line.split(" ")
                    if(len(line)>2):
                        if(line[1]=="Total"):
                            count+=1
                            # timeFTMSum+=float(line[-1])
                            # print(line[-1])
                            timeFTM+=float(line[-1])

                values.append(timeFTM)
                f.close()
        values.sort()
        values = values[1:-1]
        mean = sum(values)/len(values)
        perfTable[i+1].append(format(mean, '.2f'))

        # FTM CP
        values=[]
        for itry in range(0,12):
            timeTTK=0
            # print("./outputs_average_cp/ttk_"+dataset+para+"_try"+str(itry))
            with open("./outputs_cp/ttk_"+dataset+para+"_try"+str(itry), "r") as f:
                for line in f.readlines():
                    line=line.split(" ")
                    if(len(line)>4):
                        if(line[4]=="processed"):
                            timeTTK+=float(line[6])

                values.append(timeTTK)
                f.close()
        values.sort()
        values = values[1:-1]
        mean = sum(values)/len(values)
        perfTableCP[i+1].append(format(mean, '.2f'))

        # nonprog 
        timeNPSum=0
        values=[]
        for itry in range(0,12):
            timeNP=0
            with open("./outputs_pd//nonProg_"+dataset+para+"_try"+str(itry), "r") as f:
                for line in f.readlines():
                    line=line.split(" ")
                    if(len(line)>1):
                        line2 = line
                        if(line2[0]=='Complete'):
                            result=line2[1]
                            result=result.split("(")[1]
                            result=result.split("s")[0]
                            timeNP+=float(result)
                            values.append(timeNP)
                f.close()

        values.sort()
        values = values[1:-1]
        mean = sum(values)/len(values)
        perfTable[i+1].append(format(mean, '.2f'))

        #nonprog CP
        values=[]
        for itry in range(0,12):
            timeTTK=0
            # print("./outputs_average_cp/nonProg_"+dataset+para+"_try"+str(itry))
            with open("./outputs_cp/nonProg_"+dataset+para+"_try"+str(itry), "r") as f:
                for line in f.readlines():
                    line=line.split(" ")
                    # print(line)
                    if(len(line)>4):
                        if(line[4]=="processed"):
                            timeTTK+=float(line[6])
                            print("got ", timeTTK)

                values.append(timeTTK)
                f.close()
        values.sort()
        values = values[1:-1]
        mean = sum(values)/len(values)
        perfTableCP[i+1].append(format(mean, '.2f'))
        print(values)
        print(perfTableCP)

        # Prog
        values=[]
        for itry in range(0,12):
            time=0
            with open("./outputs_pd/prog_"+dataset+para+"_try"+str(itry),"r") as f:
                for line in f.readlines():
                    line=line.split(" ")
                    if(len(line)>1):
                        line2 = line
                        if(line2[0]=='Complete'):
                            result=line2[1]
                            result=result.split("(")[1]
                            result=result.split("s")[0]
                            time=float(result)
                            values.append(time)
                f.close()

        values.sort()
        values = values[1:-1] #removing outliers
        mean = sum(values)/len(values)
        perfTable[i+1].append(format(mean, '.2f'))

        #prog CP
        values=[]
        for itry in range(0,12):
            timeTTK=0
            with open("./outputs_cp/prog_"+dataset+para+"_try"+str(itry), "r") as f:
                for line in f.readlines():
                    line=line.split(" ")
                    if(len(line)>4):
                        if(line[4]=="processed"):
                            timeTTK+=float(line[6])

                values.append(timeTTK)
                f.close()
        values.sort()
        values = values[1:-1]
        mean = sum(values)/len(values)
        perfTableCP[i+1].append(format(mean, '.2f'))

    perfTable=np.array(perfTable)
    perfTableCP=np.array(perfTableCP)



    if(para==""):
        perfTableSeq=perfTable
        perfTableCPSeq=perfTableCP
    else:
        perfTableCPPara=perfTableCP
        perfTablePara=perfTable


table=[[""]*9]*perfTable.shape[0]
table=np.array(table, dtype='object')
table[:,0]=np.array([text[0].upper()+text[1:] for text in perfTableSeq[:,0]])
table[:,1]=perfTableCPSeq[:,1]
table[:,2]=perfTableCPSeq[:,2]
table[:,3]=perfTableCPSeq[:,3]
table[:,5]=perfTableSeq[:,1]
table[:,6]=perfTableSeq[:,2]
table[:,7]=perfTableSeq[:,3]
table[1:,4]=perfTableCPSeq[1:,2].astype(float)/perfTableCPSeq[1:,3].astype(float) 

for i in range(len(table[1:,4])):
    x=table[1:,4][i]
    table[1:,4][i]=format(x, '.2f')

table[1:,8]=perfTableSeq[1:,2].astype(float)/perfTableSeq[1:,3].astype(float) 
for i in range(len(table[1:,8])):
    x=table[1:,8][i]
    table[1:,8][i]=format(x, '.2f')


with open("./table2.tex","w") as f:
    f.write("\\begin{tabular}{|l|rrrr|rrrr|}\n")
    # f.write("%s\\\\\n" % " & ".join(str(col).title() for col in table[0]))
    f.write("\hline\n")
    f.write("    &\multicolumn{4}{c|}{Critical Points}&\multicolumn{4}{c|}{Persistence Diagram}\\\\\n")
    f.write(" Data set &TTK&NP&Prog &{\\small Speedup}& TTK & NP &  Prog & {\small Speedup}\\\\\n")
    # f.write("{Data set} & 1th & 8th & 1th & 8th& 1th & 8th\\\\\n")
    f.write("\hline\n")
    for row in table[1:]:
        f.write("%s\\\\\n" % " & ".join(col for col in row))
    f.write("\hline\n")
    f.write("\\end{tabular}\n")

    f.close()


table=[[""]*9]*perfTable.shape[0]
table=np.array(table, dtype='object')
table[:,0]=np.array([text[0].upper()+text[1:] for text in perfTablePara[:,0]])
table[:,1]=perfTableCPPara[:,1]
table[:,2]=perfTableCPPara[:,2]
table[:,3]=perfTableCPPara[:,3]
table[:,5]=perfTablePara[:,1]
table[:,6]=perfTablePara[:,2]
table[:,7]=perfTablePara[:,3]
table[1:,4]=perfTableCPPara[1:,2].astype(float)/perfTableCPPara[1:,3].astype(float) 

for i in range(len(table[1:,4])):
    x=table[1:,4][i]
    table[1:,4][i]=format(x, '.2f')

table[1:,8]=perfTablePara[1:,2].astype(float)/perfTablePara[1:,3].astype(float) 
for i in range(len(table[1:,8])):
    x=table[1:,8][i]
    table[1:,8][i]=format(x, '.2f')

with open("./table3.tex","w") as f:
    f.write("\\begin{tabular}{|l|rrrr|rrrr|}\n")
    # f.write("%s\\\\\n" % " & ".join(str(col).title() for col in table[0]))
    f.write("\hline\n")
    f.write("    &\multicolumn{4}{c|}{Critical Points}&\multicolumn{4}{c|}{Persistence Diagram}\\\\\n")
    f.write(" Data set &TTK&NP&Prog &{\\small Speedup}& TTK & NP &  Prog & {\small Speedup}\\\\\n")
    # f.write("{Data set} & 1th & 8th & 1th & 8th& 1th & 8th\\\\\n")
    f.write("\hline\n")
    for row in table[1:]:
        f.write("%s\\\\\n" % " & ".join(col for col in row))
    f.write("\hline\n")
    f.write("\\end{tabular}\n")

    f.close()
