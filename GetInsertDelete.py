# INDEL analysis method from Sean Simmons -- KK edit
# conda activate pysv


import pysam
import pandas as pd
import os
import glob




def GetInsertDel_all(bam_list,reg_df):
    regs=[[reg_df.iloc[i][0],reg_df.iloc[i][1],reg_df.iloc[i][2],reg_df.iloc[i][3]] for i in range(0,reg_df.shape[0])]
    dat=GetInsertDel_regions(bam_list,regs)
    return(dat)


##
##Gets insertion and deletions for all bams in directory, excluding those that aren't from perturbations
##Use "../Ch[1-5]" as dir
##
def GetInsertDel_directory(bams,chrom,start,end, fully_spanning=False):
    #os.system("ls -1 "+dir+"/*bam | grep -v reads | grep -v combined > temp.txt")
    #bams=open("temp.txt","r").read().strip().split("\n") ##Get list of bams
    print(bams)
    results=[GetInsertDel(bamfile,chrom,start,end,fully_spanning) for bamfile in bams]
    for i in range(0,len(bams)):
        results[i]["Name"]=["Ch"+str(i+1)+"_"+nam for nam in results[i]["CBC"]]
        results[i]["File"]=bams[i]
    dat=pd.concat(results,axis=0)
    return(dat)





##
##Runs GetInsertDel_directory for multiple regions
##
def GetInsertDel_regions(bam_list,regs, fully_spanning=False):
    results=[GetInsertDel_directory(bam_list,reg[0],reg[1],reg[2], fully_spanning) for reg in regs]
    for i in range(0,len(regs)):
        reg=regs[i]
        results[i]["Region"]=reg[0]+":"+str(reg[1])+"-"+str(reg[2])
        results[i]["Guide"]=reg[3]
    dat=pd.concat(results,axis=0)
    return(dat)

##
##Goes read by read and determines if has insertion/deletion in region of interest
##
def GetInsertDel(bamfile,chrom,start,end, fully_spanning=False):
    print(chrom+":"+str(start)+"-"+str(end))
    samfile = pysam.AlignmentFile(bamfile, "rb") ##Load sam file
    numTot=0 ##Number reads
    numInsert=0 ##Number with insertion
    numDelete=0 ##Number with deletion
    numBoth=0;
    ret=[]
    for seq in samfile.fetch(contig=chrom,start=start,end=end):
        if not seq.has_tag("CB"):
            continue;
        cbc=seq.get_tag("CB")
        [insert,delete,overlaps]=GetInsertDelete(seq,chrom,start,end, fully_spanning=fully_spanning)
        if overlaps==0:
            continue;
        numTot=numTot+1
        numInsert=numInsert+insert
        numDelete=numDelete+delete
        if insert>0 and delete>0:
            numBoth=numBoth+1
        ret.append([insert,delete,cbc])
    samfile.close()
    print(bamfile)
    print("Total "+str(numTot))
    print("Insert "+str(numInsert))
    print("Delete "+str(numDelete))
    print(" ")
    if numTot==0:
#    	ret=[[0,0,"hi"]]
#    dat=pd.DataFrame(ret)
#    dat.columns=["Insert","Delete","CBC"]
#    dat["Reads"]=1
#    dat=dat.groupby("CBC").sum().reset_index()
#    return(dat)
        return pd.DataFrame(columns=["CBC", "Insert", "Delete", "Reads"])
    dat = pd.DataFrame(ret, columns=["Insert", "Delete", "CBC"])
    dat_agg = dat.groupby("CBC").agg(
        Insert=('Insert', 'sum'),
        Delete=('Delete', 'sum'),
        Reads=('CBC', 'size')).reset_index()  # Count the number of rows (reads) in each group
    return dat_agg

##
##Counts number of deletions/insertions based on refPos
##Returns insert (0 or 1), delete (0 or 1), and overlap (0 or 1) where 1 is true, 0 false
##
def GetInsertDelete(seq,chrom,start,end, fully_spanning=False):
    if fully_spanning:
        if seq.is_unmapped or not (seq.reference_start <= start and seq.reference_end >= end):
            return [0, 0, 0]
    if chrom!=seq.reference_name:
        print("Yuck!")
        print(seq.reference_name)
        return([0,0,0])
    refPos=seq.get_reference_positions(full_length=True) ##the position each base maps to
    insert=0;
    delete=0;
    overlaps=0;
    posprev=None
    pos=None
    inIns=False
    for cur in refPos:
        posprev=pos
        pos=cur
        if pos!=None and pos>start and pos<end:
            if inIns and overlaps==1:
                insert=1;
                inIns=False
            overlaps=1;
            if posprev!=None:
                if pos>posprev+1 and pos<posprev+10:
                    delete=1
        if pos==None:
            inIns=True;
        if pos!=None:
            inIns=False
    return([insert,delete,overlaps])



# Main driver function
# 1. Finds all BAM files in the current directory
# 2. Reads target regions from "gRNA.csv"
# 3. Runs the analysis across all files and regions
# 4. Saves the consolidated results
def main():
    print("Start!")
    print("Get BAM files!")
    bam_files = glob.glob("*.bam")
    if not bam_files:
        print("Error: No .bam files found in the current directory.")
        return
    print(f"Found {len(bam_files)} BAM files to analyze.")

    print("Read CSV file to get target regions!")
    grna_data = pd.read_csv("gRNA.csv", header=None)
    grna_data.columns = ["name", "seq", "chrom", "start", "end"]
    # Format: [chrom, start, end, guide_name]
    regions = grna_data[["chrom", "start", "end", "name"]].values.tolist()
    print(f"Loaded {len(regions)} regions from gRNA.csv.")

    print("Find INDEL!")
    final_results = GetInsertDel_regions(bam_files, regions, fully_spanning=True)

    print("save!")
    output_filename = "results_fully_spanning.csv"
    final_results.to_csv(output_filename, index=False)
    print("DONE")




if __name__ == "__main__":
    main()

