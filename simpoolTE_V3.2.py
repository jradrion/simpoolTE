##Simulates both the insertion and deletion of TEs in a population of chromosomes according to the neutral frequecy distribution, then simulates pool-seq sequencing on these chromosomes

import os, sys, argparse, random, glob, numpy, copy
from operator import itemgetter
import multiprocessing as mp

def progress_bar(percent, barLen = 50):
    sys.stdout.write("\r")
    progress = ""
    for i in range(barLen):
        if i < int(barLen * percent):
            progress += "="
        else:
            progress += " "
    sys.stdout.write("[ %s ] %.2f%%" % (progress, percent * 100))
    sys.stdout.flush()

def fastaformat(seq):
    '''returns a string where every 70 bases is separated by a new line'''
    return '\n'.join(seq[i:i+70] for i in xrange(0,len(seq),70))

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]

def invertLines(samples,nChroms):
    allSamples=[]
    for x in range(nChroms):
        allSamples.append(x)
    invertedSamples=[]
    for x in allSamples:
        if not x in samples:
            invertedSamples.append(x)
    return invertedSamples

def subfinder(mylist, pattern):
        matches = []
        for i in range(len(mylist)):
            if mylist[i] == pattern[0] and mylist[i:i+len(pattern)] == pattern:
                matches.append(mylist[i:(i+len(pattern)+7)])
        return matches

def simulate(seq,teFasta,inserts,deletions,name,out2):
    '''given a sequence to simulate (in the form of a list) and the new_insertlist and new_delList for that sequence,
    write a simulated sequence'''
    for x in deletions:
        for i in xrange(x[0]-1,x[1]):
            seq[i]='%'
    count=0
    for i in xrange(len(inserts)):
        seq.insert(inserts[i][0]+count,'$')
        count+=1
    seq=remove_values_from_list(seq,'%')
    count=0
    for i in xrange(len(seq)):
        if seq[i] == '$':
            teSeq=''
            for te in teFasta:
                if inserts[count][4] in te:
                    teSeq=teFasta[te]
            tsd=''.join(seq[i-inserts[count][-1]:i])
            teSeq=teSeq+tsd
            seq[i]=teSeq
            count+=1
    joinedSeq=''.join(seq)
    with open(out2+name+'_sim.fa', 'w') as fOUT:
        fOUT.write('>'+name+'\n'+fastaformat(joinedSeq)+'\n')
    with open(out2+name+'_insert_log.txt', 'w') as fOUT:
        fOUT.write('insert_pos\tname\tlength\tfreq\tlenTSD\n')
        for x in inserts:
            fOUT.write(str(x[0])+'\t'+x[4]+'\t'+str(x[3]-x[2]+1)+'\t'+str(x[-2])+'\t'+str(x[-1])+'\n')
    with open(out2+name+'_deletion_log.txt', 'w') as fOUT:
        fOUT.write('start_pos\tend_pos\tname\tlength\tfreq\n')
        for x in deletions:
            fOUT.write(str(x[0])+'\t'+str(x[1])+'\t'+str(x[2])+'\t'+str(x[1]-x[0]+1)+'\t'+str(x[-1])+'\n')

def catPool(path,direction,out4):
    catString=''
    for i in xrange(len(path)):
        catString+=path[i]+' '
    os.system('cat '+catString+'> '+ out4 + 'simPool_' + direction + '.fq')
    print 'Pool',direction,'completed!'

def rehead(path,identifier):
    with open(path,'r') as fIN, open(path.replace('.fq','_rehead.fq'), 'w') as fOUT:
        for line in fIN:
            if '@' in line:
                arr=line.split('/')
                line=arr[0]+'_'+identifier+'/'+arr[1]
                fOUT.write(line)
            else:
                fOUT.write(line)
    print identifier,'rehead complete!'

def isNested(chrm,start,stop,ID,annotation,L):
    '''return 1 if focal TE is nested/overlaps within another reference TE
    TEs are considered overlapping if they are < averege read length from another TE'''
    for x in annotation:
        if ID != x[3] and chrm==x[0] and x[1]-L < start and stop < x[2]+L:
            return 1
        if ID != x[3] and chrm==x[0] and x[1]-L < start and start < x[2]+L:
            return 1
        if ID != x[3] and chrm==x[0] and x[1]-L < stop and stop < x[2]+L:
            return 1
    return 0

def containsNested(chrm,start,stop,ID,annotation,L):
    '''return 1 if focal TE contains another reference TE nested within itself'''
    for x in annotation:
        if ID != x[3] and chrm==x[0] and start-L < x[1] and x[2] < stop+L:
            return 1
    return 0

def generateTSD(lam):
    '''generate a random tsd length based on a poisson distribution around lambda'''
    d=numpy.random.poisson(lam,1000)
    return d[random.randint(0,999)]

def removeBedPos(bed, genome):
    extractedSeqs={}
    for ch in genome:
        for x in bed:
            if ch == x[0]:
                extractedSeqs[x[3]]=""
                for i in xrange(x[1]-1,x[2]):
                    extractedSeqs[x[3]]+=genome[ch][i]
    return extractedSeqs

def assign_task(siteID, task_q, nProcs):
    c,i,nth_job=0,0,1
    while (i+1)*nProcs <= len(siteID):
        i+=1
    nP1=nProcs-(len(siteID)%nProcs)
    for j in range(nP1):
        task_q.put((siteID[c:c+i], nth_job))
        nth_job += 1
        c=c+i
    for j in range(nProcs-nP1):
        task_q.put((siteID[c:c+i+1], nth_job))
        nth_job += 1
        c=c+i+1

def create_proc1(nProcs, task_q, params):
    for _ in range(nProcs):
        p = mp.Process(target=worker1, args=(task_q, params))
        p.daemon = True
        p.start()

def create_proc2(nProcs, task_q, params):
    for _ in range(nProcs):
        p = mp.Process(target=worker2, args=(task_q, params))
        p.daemon = True
        p.start()

def create_proc3(nProcs, task_q, params):
    for _ in range(nProcs):
        p = mp.Process(target=worker3, args=(task_q, params))
        p.daemon = True
        p.start()

def worker1(task_q, params):
    while True:
        try:
            groups, nth_job = task_q.get()
            #unpack parameters
            pirsPATH,refChrom,out1 = params
            for i in groups:
                cmd=pirsPATH +' diploid -i '+ refChrom + ' -d 0.0 -v 0.0 -c 0 -o ' + out1 + 'sample.' + str(i+1)
                os.system(cmd)
        finally:
            task_q.task_done()

def worker2(task_q, params):
    while True:
        try:
            groups, nth_job = task_q.get()
            #unpack parameters
            chromosomes,teFasta,new_insertList,new_delList,out2 = params
            for i in groups:
                simulate(chromosomes[i],teFasta,new_insertList[i],new_delList[i],'sample.' + str(i+1),out2)
                print 'Sample',i+1,'simulated!'
        finally:
            task_q.task_done()

def worker3(task_q, params):
    while True:
        try:
            paths, nth_job = task_q.get()
            #unpack parameters
            pirsPATH,rLen,cov,insz,out3 = params
            for path in paths:
                cmd=pirsPATH +' simulate -i '+ path + ' -l '+str(rLen)+' -x '+str(cov)+' -m '+str(insz)+' -c 0 -o ' + out3 + path.split('/')[-1].split('.fa')[0]
                os.system(cmd)
        finally:
            task_q.task_done()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory',default=-1)
    parser.add_argument('-pirs',dest='pirsPATH',help='pirs path')
    parser.add_argument('-seqtk',dest='seqtkPATH',help='seqtk path')
    parser.add_argument('-c',dest='refChrom',help='chromosome to simulate fasta format')
    parser.add_argument('-b',dest='bed',help='te annotation bed file')
    parser.add_argument("-ex",dest="exclude",help="newline separated list of TEs to exclude from simulation(name must match that from bed file)",default="none")
    parser.add_argument('-mnlen',dest='minlength',help='minimum length of TEs to insert and delete', type=int, default=500)
    parser.add_argument('-mxlen',dest='maxlength',help='maximum length of TEs to insert and delete', type=int, default=10000)
    parser.add_argument('-nchr',dest='numChroms',help='number of chromosomes in the population', type=int, default=10)
    parser.add_argument('-nte',dest='numTEs',help='number of insertions and deletion to simulate in lowest frequency class', type=int, default=50)
    parser.add_argument('-r',dest='randseed',help='seed for random number generator', type=int, default=12345)
    parser.add_argument('-x',dest='coverage',help='coverage to simulate',type=int, default='50')
    parser.add_argument('-rlen',dest='readlen',help='read length to simulate',type=int, default='100')
    parser.add_argument('-insz',dest='insz',help='insert size to simulate',type=int,  default='200')
    parser.add_argument('-t',dest='threads',help='number of threads',  type=int, default=1)
    args = parser.parse_args()

    # identify current working directory
    if args.wd == -1:
        cwd=os.getcwd()
    else:
        cwd=os.path.realpath(args.wd)

    #Import options
    mnLen=args.minlength
    mxLen=args.maxlength
    nChroms=args.numChroms
    nTEs=args.numTEs
    randSeed=args.randseed
    cov=args.coverage
    rLen=args.readlen
    insz=args.insz
    nProcs=args.threads

    #Initialize random number generator
    random.seed(randSeed)

    #Read bed file
    bedEntries = []
    with open(args.bed, 'r') as infile:
        for line in infile:
            arr=line.split()
            bedEntries.append([arr[0],int(arr[1]),int(arr[2]),arr[3],arr[5]])

    #Allocate indel counts and frequencies for each bin
    dist_fq,dist_counts=[],[]
    for x in range(nChroms-1):
        dist_fq.append(round(1/(float(x)+1),3))
    for i in xrange(len(dist_fq)):
        dist_counts.append(int(round(dist_fq[i]*nTEs,0)))
    total_inserts=sum(dist_counts)

    #Read teChromosomes to exclude
    if args.exclude == "none":
        excludeList=0
    else:
        excludeList=[]
        with open(args.exclude, "r") as fIN:
            for line in fIN:
                excludeList.append(line.split()[0])

    #Check to make sure there are enough refernce tes to delete
    with open(args.refChrom, "r") as fIN:
        for line in fIN:
            if line.startswith(">"):
                chromosome=line.split()[0][1:]
                break
    count=0
    for x in bedEntries:
        if x[0] == chromosome and mxLen > (x[2] - x[1]) >= mnLen and isNested(x[0],x[1],x[2],x[3],bedEntries,rLen) == 0 and containsNested(x[0],x[1],x[2],x[3],bedEntries,rLen) == 0:
            if excludeList==0:
                count +=1
            else:
                if x[3] not in excludeList:
                    count+=1
    if total_inserts > count:
        print 'WARNING: Too few reference TEs available to delete given current parameters!'
        print 'Attempting to simulate',total_inserts,'TE deletions with only',count,'reference TEs available to delete!'
        print 'Please reduce either the number of chromosomes/TEs or change the min and max length thresholds!'
        sys.exit()
    else:
        print "Simulating the insertion and deletion of %s unique TEs..." %(total_inserts)

    #Create new directories
    os.system('mkdir '+cwd+'pirs_snps')
    os.system('mkdir '+cwd+'teSim_simData')
    os.system('mkdir '+cwd+'pirs_reads')
    os.system('mkdir '+cwd+'teSim_pooledReads')
    out1=cwd+'pirs_snps/'
    out2=cwd+'teSim_simData/'
    out3=cwd+'pirs_reads/'
    out4=cwd+'teSim_pooledReads/'

    #Read chromosome to simulate
    print "Reading chromosome to simulate..."
    validBase="ACGT"
    invalidIndex,ct=[],0
    rawChr,acgtN = "",[0,0,0,0,0]
    with open(args.refChrom, "r") as fIN:
        for line in fIN:
            if line.startswith(">"):
                chromosome=line.split()[0][1:]
            else:
                seq=line.replace("\n","").upper()
                for s in seq:
                    if s in validBase:
                        acgtN[validBase.index(s)]+=1
                        ct+=1
                    else:
                        acgtN[-1]+=1
                        invalidIndex.append(ct)
                        ct+=1
                rawChr+=seq

    #Draw ACGT from their density distributions
    print "Drawing ACGT characters with from their density distributions..."
    newBase=[]
    for i in range(len(validBase)):
        n=int((acgtN[i]/float(sum(acgtN[:4])))*acgtN[-1]+1)
        print "Drawing %s's..." %(validBase[i])
        for j in range(n):
            progress_bar(j/float(n))
            newBase.append(validBase[i])
        print "\n"
    random.shuffle(newBase)

    #Replace non-ACGT characters with ACGT's
    print "Replacing non-ACGT characters with ACGT's..."
    ct=0
    rawChrList=list(rawChr)
    for i in range(len(invalidIndex)):
        progress_bar(i/float(len(invalidIndex)))
        rawChrList[invalidIndex[i]]=newBase[i]
    print "\n"
    rawChr="".join(rawChrList)

    #Extract annotated bed sequences to te dictionary with fasta seqs
    refChrom={}
    refChrom[chromosome]=rawChr
    teFasta=removeBedPos(bedEntries,refChrom)

    #Run multiprocess 1
    print "Simulating SNPs..."
    task_q = mp.JoinableQueue()
    params=[args.pirsPATH,args.refChrom,out1]
    create_proc1(nProcs, task_q, params)
    assign_task(range(nChroms), task_q, nProcs)
    try:
        task_q.join()
    except KeyboardInterrupt:
        print "KeyboardInterrupt"
        sys.exit(0)
    else:
        print "Finished!"

    #Read the sequences from the input directory
    paths = []
    for file in glob.glob(os.path.join(out1, '*.fa')):
        ct=int(file.split(".")[-3])
        paths.append([ct,file])

    #Read chromosome into list
    print "\nReading and converting simulated samples..."
    chromosomes=[]
    for i in range(len(paths)):
        progress_bar(i/float(len(paths)-1))
        raw=""
        with open(sorted(paths)[i][1], "r") as fIN:
            for line in fIN:
                if not line.startswith(">"):
                    seq=line.replace("\n","").upper()
                    raw+=seq
        rawList=list(raw)
        chromosomes.append(rawList)
    print "\n"

    #Randomly select TEs to insert
    print "Selecting TEs to insert..."
    inPool,icheck,ct=[],[],0
    while len(inPool) < total_inserts:
        progress_bar(ct/(float(total_inserts)-1))
        i=random.randint(0,len(bedEntries)-1)
        if not i in icheck:
            if bedEntries[i][0] == chromosome and args.minlength <= (bedEntries[i][2] - bedEntries[i][1]) <= mxLen:
                if excludeList == 0:
                    if containsNested(bedEntries[i][0],bedEntries[i][1],bedEntries[i][2],bedEntries[i][3],bedEntries,rLen) == 0:
                        inPool.append(bedEntries[i])
                        icheck.append(i)
                        ct+=1
                else:
                    if not bedEntries[i][3] in excludeList:
                        if containsNested(bedEntries[i][0],bedEntries[i][1],bedEntries[i][2],bedEntries[i][3],bedEntries,rLen) == 0:
                            inPool.append(bedEntries[i])
                            icheck.append(i)
                            ct+=1
    print "\n"

    #Randomly select TEs to delete
    print "Selecting TEs to delete..."
    random.seed(randSeed+1)
    delPool,icheck,ct=[],[],0
    while len(delPool) < total_inserts:
        progress_bar(ct/(float(total_inserts)-1))
        i=random.randint(0,len(bedEntries)-1)
        if not i in icheck:
            if bedEntries[i][0] == chromosome and args.minlength <= (bedEntries[i][2] - bedEntries[i][1]) <= mxLen:
                if excludeList == 0:
                    if isNested(bedEntries[i][0],bedEntries[i][1],bedEntries[i][2],bedEntries[i][3],bedEntries,rLen) == 0 and containsNested(bedEntries[i][0],bedEntries[i][1],bedEntries[i][2],bedEntries[i][3],bedEntries,rLen) == 0:
                        delPool.append(bedEntries[i])
                        icheck.append(i)
                        ct+=1
                else:
                    if not bedEntries[i][3] in excludeList:
                        if isNested(bedEntries[i][0],bedEntries[i][1],bedEntries[i][2],bedEntries[i][3],bedEntries,rLen) == 0 and containsNested(bedEntries[i][0],bedEntries[i][1],bedEntries[i][2],bedEntries[i][3],bedEntries,rLen) == 0:
                            delPool.append(bedEntries[i])
                            icheck.append(i)
                            ct+=1
    print "\n"
    tmp=copy.deepcopy(delPool)
    delPool=tmp

    #Randomly select in which lines to insert and delete
    index,insertLine,counter=[],[],0
    for i in xrange(len(dist_counts)):
        temp=[]
        for j in xrange(dist_counts[i]):
            temp.append(counter)
            counter+=1
        index.append(temp)
    for i in xrange(len(inPool)):
        for j in xrange(len(index)):
            for l in xrange(len(index[j])):
                if i == index[j][l]:
                    insertLine.append(random.sample(range(nChroms),j+1))
    ct=0
    for x in inPool:
        x.append(insertLine[ct])
        ct+=1
    ct=0
    for x in delPool:
        x.append(invertLines(insertLine[ct],nChroms))
        ct+=1

    #Randomly select sites to insert that do not overlap with annotated TEs
    print "Identifying insertion sites..."
    lens=[]
    for x in chromosomes:
        lens.append(len(x))
    insert_pos,insert_check,ct=[],[],0
    while len(insert_pos) < len(inPool):
        progress_bar(ct/(float(len(inPool))-1))
        x=random.randint(rLen,max(lens)-rLen)
        for i in xrange(len(inPool)):
            if isNested(chromosome,x,x,"!novel!",bedEntries,rLen) == 0 and not x in insert_check:
                insert_check.append(x)
                insert_pos.append(x)
                ct+=1
    print "\n"
    for i in xrange(len(inPool)):
        fq=len(inPool[i][-1])/float(nChroms)
        inPool[i].append(insert_pos[i])
        inPool[i].append(fq)
        inPool[i].append(generateTSD(5))
    for i in range(len(delPool)):
        fq=1-(len(delPool[i][-1])/float(nChroms))
        delPool[i].append(fq)

    #Write the list of simulated insertion and deletions
    #inPool=[chrm,start,stop,name,strand,samples,insertion position, fq, tsd]
    #delPool=[chrm,start,stop,name,strand,samples,fq]
    with open(out2+'simulated_insertion_list.txt', 'w') as fOUT:
        fOUT.write('chr\tinsert_site\tname\tlength\tfreq\tlenTSD\tsamples\n')
        for x in sorted(inPool, key=itemgetter(6)):
            fOUT.write(x[0]+'\t'+str(x[6])+'\t'+str(x[3])+'\t'+str(x[2]-x[1]+1)+'\t'+str(x[7])+'\t'+str(x[8])+'\t'+str(x[5])+'\n')
    with open(out2+'simulated_deletion_list.txt', 'w') as fOUT:
        fOUT.write('chr\tstart\tend\tname\tlength\tfreq\tsamples\n')
        for x in sorted(delPool, key=itemgetter(1)):
            fOUT.write(x[0]+'\t'+str(x[1])+'\t'+str(x[2])+'\t'+str(x[3])+'\t'+str(x[2]-x[1]+1)+'\t'+str(x[6])+'\t'+str(x[5])+'\n')

    #Create a list of tes that will be inserted and deleted
    insert_list,del_list=[],[]
    for i in xrange(nChroms):
        temp_in,temp_del=[],[]
        for x in sorted(inPool, key=itemgetter(6)):
            if i in x[5]:
                temp_in.append(x)
        insert_list.append(temp_in)
        for x in sorted(delPool, key=itemgetter(1)):
            if i in x[5]:
                temp_del.append(x)
        del_list.append(temp_del)
    new_insertList,new_delList=[],[]
    for x in insert_list:
        temp=[]
        for i in xrange(len(x)):
            temp.append([x[i][6],x[i][0],x[i][1],x[i][2],x[i][3],x[i][7],x[i][8]])
        new_insertList.append(temp)
    for x in del_list:
        temp=[]
        for i in xrange(len(x)):
            temp.append([x[i][1],x[i][2],x[i][3],x[i][4],x[i][6]])
        new_delList.append(temp)

    #Run multiprocess 2
    print "Simulating TE insertions and deletions..."
    task_q = mp.JoinableQueue()
    params=[chromosomes,teFasta,new_insertList,new_delList,out2]
    create_proc2(nProcs, task_q, params)
    assign_task(range(nChroms), task_q, nProcs)
    try:
        task_q.join()
    except KeyboardInterrupt:
        print "KeyboardInterrupt"
        sys.exit(0)
    else:
        print "Finished!"

    #Simulate sequencing
    print 'Simulating paired-end sequencing...'
    paths = []
    for file in glob.glob(os.path.join(out2, '*.fa')):
        paths.append([file.split(".")[-3],file])
    sortedPaths=[]
    for x in sorted(paths):
        sortedPaths.append(x[1])

    #Run multiprocess 3
    task_q = mp.JoinableQueue()
    params=[args.pirsPATH,rLen,cov,insz,out3]
    create_proc3(nProcs, task_q, params)
    assign_task(sortedPaths, task_q, nProcs)
    try:
        task_q.join()
    except KeyboardInterrupt:
        print "KeyboardInterrupt"
        sys.exit(0)
    else:
        print "Finished!"

    #Randomly sample sequences to create pool
    print 'Sampling reads to pool...'
    paths = []
    for file in glob.glob(os.path.join(out3, '*.fq')):
        paths.append(file)
    paths.sort()
    reads=1/(len(paths)*0.5)
    seeds=[]
    counter=0
    with open(out4+'pooling_seeds.txt', 'w') as fOUT:
        for i in range(len(paths)):
            if '_1.fq' in paths[i]:
                seed=random.randint(0,1000)
                seeds.append(seed)
                os.system(args.seqtkPATH +' sample -s '+str(seed)+' '+ paths[i] + ' ' + str(reads) + ' > ' + out3 + paths[i].split('/')[-1].replace('.fq','_pool_sampled.fq'))
                fOUT.write(args.seqtkPATH +' sample -s '+str(seed)+' '+ paths[i] + ' ' + str(reads) + ' > ' + out3 + paths[i].split('/')[-1].replace('.fq','_pool_sampled.fq')+'\n')
        for i in range(len(paths)):
            if '_2.fq' in paths[i]:
                seed=seeds[counter]
                os.system(args.seqtkPATH +' sample -s '+str(seed)+' '+ paths[i] + ' ' + str(reads) + ' > ' + out3 + paths[i].split('/')[-1].replace('.fq','_pool_sampled.fq'))
                fOUT.write(args.seqtkPATH +' sample -s '+str(seed)+' '+ paths[i] + ' ' + str(reads) + ' > ' + out3 + paths[i].split('/')[-1].replace('.fq','_pool_sampled.fq')+'\n')
                counter+=1
    paths=[]
    for file in glob.glob(os.path.join(out3, '*_pool_sampled.fq')):
        paths.append(file)
    paths.sort()
    for file in paths:
        rehead(file,file.split('/')[-1].split('_')[0])
    print 'Generating pool...'
    paths=[]
    for file in glob.glob(os.path.join(out3, '*pool_sampled_rehead.fq')):
        paths.append(file)
    paths.sort()
    forward,reverse=[],[]
    for x in paths:
        if int(x.split('_')[-4]) % 2 == 0:
            reverse.append(x)
        else:
            forward.append(x)
    catPool(forward,'1',out4)
    catPool(reverse,'2',out4)

    #Write new chromosome for reference
    with open(os.path.join(out2,"refChrom_ACGT_replaced.fa"), "w") as fOUT:
        fOUT.write(">%s\n" %(chromosome))
        fOUT.write(fastaformat(rawChr)+"\n")

    print "teSim_V3.1 finished!"


if __name__ == "__main__":
    main()

