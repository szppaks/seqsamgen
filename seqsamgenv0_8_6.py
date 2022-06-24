#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 21:09:13 2022

@author: Péter Szutor, Attila Horváth
v0.8.6
"""
import numpy.random as nprand
import numpy as np
from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import seqlogo
import sys
import os
import pandas
import tkinter as tk
from tkinter import ttk
from io import StringIO
from tkinter import messagebox
import configparser
import string
nuks=['A','C','G','T']
nuks2=['G','A','T','C']

def homerjavito(fhomer):
    sor=fhomer.readline()
    ski=''
    while sor:
        if sor[0]!='>' and len(sor)>0:
            szamok=sor.replace('\t',' ').replace('   ',' ').replace('  ',' ')
            t=np.fromstring(szamok,sep=' ')
            if np.sum(t)!=1:
                diff=round((np.sum(t)-1)/len(t),5)
                t=np.round(t-diff,5)
                s='\t'.join([str(ti) for ti in t])
                ski=ski+s+'\n'
            else:
                ski=ski+sor
        else:
            ski=ski+sor
        sor=fhomer.readline()
    else:
        ski=ski+sor+'\n'
    return StringIO(ski)


def logprint(*args):
    if verbosemode:
        print(*args)

def onemotifgen(pmotif,ptip):
    try:
        ff=np.array([pmotif.counts['A',:],pmotif.counts['C',:],pmotif.counts['G',:],pmotif.counts['T',:]])
        if ptip=='pfm-four-columns':
            appm=pandas.DataFrame.from_dict(pmotif.counts)
            nuksjo=nuks2
        elif ptip=='jaspar':
            appm=seqlogo.Pfm(ff)
            appm=(seqlogo.pfm2ppm(smf))
            nuksjo=nuks
        s=''
        for i in appm.to_numpy():
            s=s+(nprand.choice(nuksjo,p=i))
    except Exception as ex:
        print(str(ex))
    return s

def motifvalidate(pmotifname):
    try:
        ki=[]
        motifext=os.path.splitext(pmotifname)[1]
        if motifext=='.jaspar':
            motip='jaspar'
        elif motifext=='.motif':
            motip='pfm-four-columns'
        else:
            print('Unknown motif format. (known formats:jaspar pfm, homer ppm)')
            return ki
        fh = open(pmotifname)
        if True:
            for m in motifs.parse(fh, motip):
                ki.append(onemotifgen(m,motip))
        retval=True
    except Exception as ex:
        print(pmotifname,' is not valid file' )
        for  onem in ki:
            print (onem.count)
        retval=False
    return retval

def motifgen(pjasparfile,pjasparname=''):
    try:
        ki=[]
        motifext=os.path.splitext(pjasparfile)[1]
        if motifext=='.jaspar':
            motip='jaspar'
            fh = open(pjasparfile)
        elif motifext=='.motif':
            motip='pfm-four-columns'
            fh = homerjavito(open(pjasparfile))
        else:
            print('Unknown motif format. (known formats:jaspar pfm, homer ppm)')
            return ki
        if pjasparname=='(ALL)' or  pjasparname=='':  #mindet hozzáaadja
            for m in motifs.parse(fh, motip, strict=False):
                ki.append(onemotifgen(m,motip))
        else:
            ml=[]
            for m in motifs.parse(fh, "jaspar"):
                if m.name==pjasparname:
                    ml.append(m)
            if len(ml)>0:
                ki.append(onemotifgen(ml[0],motip))
            else:
                print('Cannot find given record in JASPAR file')
    except Exception as ex:
        print(str(ex))
    return ki



def basednagen(pwordlen):
    s=[]
    #nuksteszt=['*',"*",'*','*']
    for i in range(pwordlen):
        s.append(nprand.choice(nuks).lower())
    return s

def motiflistlen(motiflist):
    clen=0
    for moti in motiflist:
        clen=clen+len(moti)
    return clen


def motifposgen(motifs,maxlen,gaps):
    ki=[]
    chossz=0
    hosszak=[]
    for moti in motifs:
        chossz=chossz+len(moti)
        hosszak.append(len(moti))
    hosszak=np.array(hosszak,dtype='int32')
    minkell=chossz+np.sum(gaps)
    if minkell>maxlen:
        print('Seq len too small.')
    else:
        marad=maxlen
        vegsohely=0
        for i in range(len(hosszak)):
            hely=marad-(np.sum(hosszak[:])+np.sum(gaps))
            sorsol=random.randint(0,hely)
            gaplen=0
            if len(gaps)>0 and i>0:
                gaplen=gaps[i-1]
                sorsol=gaplen
            marad=marad-(sorsol+hosszak[0])
            ki.append(vegsohely+sorsol)
            vegsohely=vegsohely+sorsol+hosszak[0]
            if len(hosszak)>1:
                hosszak=hosszak[1:]
    return ki

def gapgen(pgaps,glen,maxplace): #gap generation from possible gap values. The sum of gaps<maxplace
    gaps=pgaps[pgaps<=maxplace]
    if len(gaps)>0 and np.amin(gaps)*glen<=maxplace:
        megvan=False
        counts=0
        while not megvan and counts<10000:
            ki=[]
            for i in range(0,glen):
                ki.append(np.random.choice(gaps))
            ki=np.array(ki,dtype='int32')
            if np.sum(ki)<=maxplace:
                megvan=True
            counts+=1
    elif len(gaps)>0 and np.amin(gaps)*glen==maxplace: # pont annyi hely van, hogy csak aminimuok férnek be, akkor feltöltjük a minimummal
        ki=np.zeros((glen),dtype='int32')+np.amin(gaps)
    else:
        ki=np.array([],dtype='int32')
    return ki

def badgapgen(gaps,glen,maxgaplen):  #gap generation for bad gaps
    mindgap=np.arange(0,maxgaplen+1)
    if len(gaps)>0:
        badgaps=np.setdiff1d(mindgap,gaps)
    else:
        badgaps=mindgap
    ki=gapgen(badgaps,glen,maxgaplen)
    return ki



def seqgen(seqnum,pjasparfile,pjasparname,maxlen,pgaps,ordered='N',fileprefix=''):
    #try:
    if True:
        #create output directory
        if fileprefix is None or fileprefix=='':
            fileprefix=''.join(random.choices(string.ascii_uppercase + string.digits, k=8))

        isExist = os.path.exists(fileprefix)
        if not isExist:
            os.makedirs(fileprefix)

        gaps=np.fromstring(pgaps,dtype='int32',sep=',')
        #positive set
        kiseqs=[]
        mfs=motifgen(pjasparfile,pjasparname)
        t=[len(x) for x in mfs]
        t2=np.array(t)
        if len(gaps)>0:
            gapmin=np.amin(gaps)
        else:
            gapmin=0
        if maxlen<(np.sum(t2)+(gapmin*(len(mfs)-1))):
            print('Maximum sequence length is too small!')
            return
        for i in range(seqnum):
            mfs=motifgen(pjasparfile,pjasparname)
            if len(mfs)>0:
                s=(basednagen(maxlen)) #base random seq
                if pgaps!='':
                    goodgaps=gapgen(gaps,len(mfs)-1,maxlen-motiflistlen(mfs)) #gap generating
                else:
                    goodgaps=badgapgen(gaps,len(mfs)-1,maxlen-motiflistlen(mfs)) #gap generating
                pmfs=mfs
                if ordered=='N':
                    random.shuffle(pmfs) #if random, the shuffle
                posok=motifposgen(pmfs,maxlen,goodgaps)
                logprint('-'*maxlen)
                for j in range(len(posok)):
                    logprint('generated seq:',pmfs[j])
                    logprint('position:',posok[j])
                    es=s[posok[j]:(posok[j]+len(pmfs[j]))]=list(pmfs[j])
                seqstr=("".join(s))
                logprint(seqstr)
                s1=Seq(seqstr.upper())
                sr=SeqRecord(s1,id='pos_'+str(i),description='')
                kiseqs.append(sr)
        SeqIO.write(kiseqs, fileprefix+"/positives.fa", "fasta")

        #negative set
        #cases
        #0: random, 1: N-1 motif, 2: swapped, 3:different gap
        logprint('*'*40)
        logprint('negative set')
        logprint('*'*40)
        kiseqs=[[],[],[],[]]
        errornames=['Randomseq','n-1','shuffle','diffgap']
        proba=motifgen(pjasparfile,pjasparname)
        chossz=0
        for moti in proba:
            chossz=chossz+len(moti)
        badgaps=gapgen(gaps,len(mfs)-1,maxlen-(chossz)) #bad gap generating
        if len(proba)==1:
            errortype=lambda a : 0
            e_types=[0]
        elif ordered=='N' and (len(gaps)==0 or (len(gaps)==1 and gaps[0]==0)):
            errortype=lambda a : random.randint(0,1)
            e_types=[0,1]
        elif ordered=='Y' and (len(gaps)==0 or (len(gaps)==1 and gaps[0]==0)):
            errortype=lambda a : random.randint(0,2)
            e_types=[0,1,2]
        elif ordered=='N' and len(gaps)>0:
            errortype=lambda a : random.choice([0,1,3])
            e_types=[0,1,3]
        elif ordered=='Y' and len(gaps)>0:
            errortype=lambda a : random.choice([1,2,3])
            e_types=[0,1,2,3]
        for i in range(seqnum):
#            mitront=errortype(0)
            for mitront in e_types:
                s=(basednagen(maxlen))
                logprint('-'*maxlen)
                if mitront>0:
                    mfs=motifgen(pjasparfile,pjasparname)
                    pmfs=mfs
                    if mitront==2:
                        errname=('shuffle')
                        if len(pmfs)>2:
                            random.shuffle(pmfs)
                        elif len(pmfs)==2:
                            csere=pmfs[0]
                            pmfs[0]=pmfs[1]
                            pmfs[1]=csere
                    elif mitront==1:
                        errname=('N-1 motif')
                        ranselpos=random.randint(0,len(pmfs)-1)  #motif random erase
                        del pmfs[ranselpos]
                        ranselpos=random.randint(0,len(pmfs)-1) #gap random erase
                        if len(gaps)>1:
                            badgaps=np.delete(badgaps,ranselpos)
                    elif mitront==3:
                        errname=('Diff.gap')
                        badgaps=badgapgen(gaps,len(mfs)-1,maxlen-motiflistlen(mfs)-1) #gap generating
                    posok=motifposgen(pmfs,maxlen,badgaps)
                    logprint('-'*maxlen)
                    for j in range(len(posok)):
                        logprint('generated seq:',pmfs[j])
                        logprint('position:',posok[j])
                        es=s[posok[j]:(posok[j]+len(pmfs[j]))]=list(pmfs[j])
                else:
                    errname='Randomseq'
                logprint(errname)
                seqstr=("".join(s))
                logprint(seqstr)
                s1=Seq(seqstr.upper())
                sr=SeqRecord(s1,id='neg_'+str(i),description=errname)
                kiseqs[mitront].append(sr)

        for i in range(0,len(kiseqs)):
            if len(kiseqs[i])>0:
                SeqIO.write(kiseqs[i], fileprefix+'/'+errornames[i]+".fa", "fasta")
    # except Exception as ex:
    #     print(str(ex))
        print('Writing .fa files to '+fileprefix+' directory finished.')
        flog=open(fileprefix+'/'+fileprefix+'.log','w')
        flog.write('seqnum :%s'%str(seqnum))
        flog.write(' ; motiffilename :%s'%str(pjasparfile))
        flog.write(' ; motifname :%s'%str(pjasparname))
        flog.write(' ; gaps :%s'%str(pgaps))
        flog.write(' ; sequence length :%s'%str(maxlen))
        flog.write(' ; ordered :%s'%str(maxlen))
        flog.close()

    return

def guistart():
    root = tk.Tk()
    root.title("Sequence sample generator")
    lbl1 = tk.Label(root, text = "Number of samples (integer)")
    lbl1.grid()
    txtseqnum = tk.Entry(root, width=10)
    txtseqnum.grid(column=0)
    lbl2 = tk.Label(root, text = "Homer motif or Jaspar PFM file (filename with relative or absolute path)")
    lbl2.grid()
    txtmofile = tk.Entry(root, width=100)
    txtmofile.grid(column=0)
    lbl3 = tk.Label(root, text = "Spec motif name or empty for use all motif ")
    lbl3.grid()
    txtmoname = tk.Entry(root, width=100)
    txtmoname.grid(column=0)
    lbl4 = tk.Label(root, text = "Seq length (integer)")
    lbl4.grid()
    txtseqlen = tk.Entry(root, width=100)
    txtseqlen.grid(column=0)
    lbl5 = tk.Label(root, text = "GAP length (integer)")
    lbl5.grid()
    txtgaps = tk.Entry(root, width=100)
    txtgaps.grid(column=0)
    lbl6 = tk.Label(root, text = "Ordered")
    lbl6.grid()
    varord = tk.StringVar()
    chbord=ttk.Checkbutton(root,
                text='Order',
                variable=varord,
                onvalue='Y',
                offvalue='N')
    chbord.grid(column=0)
    lbl8 = tk.Label(root, text = "Output directory name. If not set, it will be generated automatically")
    lbl8.grid()
    txtfilepref = tk.Entry(root, width=100)
    txtfilepref.grid(column=0)
    lbl7 = tk.Label(root, text = "You can run in command mode: python seqsamgen.py [configfile.conf]")
    lbl7.grid(column=0)
    txcfg=tk.Text(bg = "light gray")
    txcfg.insert(tk.END,conftmpl)
    txcfg.grid(column=0)
    def clicked():
        try:
            seqnum=int(txtseqnum.get())
            mofile=txtmofile.get()
            moname=txtmoname.get()
            seqlen=int(txtseqlen.get())
            gaps=txtgaps.get()
            if gaps=='':
                gaps='0'
            ordering=varord.get()
            pfilepref=txtfilepref.get()
            seqgen(seqnum,mofile,moname,seqlen,gaps,ordering,pfilepref)
            messagebox.showinfo("Generating finished", "Writing to ./positives.faa and ./negatives.faa finished.")
        except:
            messagebox.showerror("Error", "Bad parameter values or missed motif file. Run in command mode for more info.")

    btn = tk.Button(root, text = "Generate seqs" ,
                 fg = "red", command=clicked)
    btn.grid(column=1,row=1)
    root.mainloop()
    return


conftmpl='''
#seqsamgen config file sample
[sequence generating]
# the number of the samples
seqnum=1000
#MOTIF file name (JASPAR PFM or HOMER PPM format)
motiffilename=example.motif
#If you want to use only one motif from file, vou can set the name.
#motifname=
#Gap sizes separated with comma.
gaps=1,3,5
#The length of the generated sequences
sequence length=120
#Motifs are ordered or not Y/N
ordered=Y

[other]
#Verbose mode: set Y if qou want to print the generated motifs
verbose mode=N
#output directory name. If not set, it will be generated automatically
output name=seqoutdir
'''

helptext='''

Sequence generator for AI training
required modules :biopython, seqlogo, pandas,numpy
It has two usage modes: command line mode with paramter file or gui mode
command line mode: python seqsamgen.py  [controlfile.cfg]
GUI mode:          python seqsamgen.py
You can set: motif file, order, gap, seq length, number of sequenties, output directory name (if not set, it generated randomly)
Known motif files: JASPAR PFM, HOMER PPM
Output format: FASTA
The control (negative) set generated automically.

Rules of generating
Motif number  |   ordered    | Gaps       |            Control set
1             |     -        |    -       |  Random
N             |     N        |    N       |  Random, N-1 motif
N             |     Y        |    N       |  Random, N-1 motif, swapped
N             |     N        |    Y       |  Random, N-1 motif, different gap
N             |     Y        |    Y       |  N-1 motif, different gap, swapped
Different gap: if you give gaps: 1,3,5 the different gap :2,4,6.. maxpossiblegaplen



'''

verbosemode=True
if len(sys.argv)<2:
    # seqgen(int(sys.argv[1]),sys.argv[2],sys.argv[3],int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]))
    # if motifvalidate("./motif2.motif"):
    #     seqgen(20,"./motif2.motif",'(ALL)',23,'','Y')
    guistart()
else:
    def getconfig(pcfg,pcname):
        try:
            ki=pcfg[pcname]
            if ki==None:
                ki=''
        except:
            ki=''
        return ki
    try:
        sgcfg = configparser.ConfigParser()
        sgcfg.read(sys.argv[1])
        c_seqsets=sgcfg['sequence generating']
        c_prgset=sgcfg['other']
        c_seqnum=int(getconfig(c_seqsets,'seqnum'))
        c_mfile=getconfig(c_seqsets,'motiffilename')
        c_mname=getconfig(c_seqsets,'motifname')
        c_gapsm=getconfig(c_seqsets,'gaps')

        c_ordered=getconfig(c_seqsets,'ordered')
        c_seqlen=int(getconfig(c_seqsets,'sequence length'))
        c_verbose=getconfig(c_prgset,'verbose mode')
        c_outputname=getconfig(c_prgset,'output name')
    except Exception as ex:
        print('Error in the config file')
        print(str(ex))
        exit()
    if motifvalidate(c_mfile):
        if c_verbose=='N':
            verbosemode=False
        seqgen(c_seqnum,c_mfile,c_mname,c_seqlen,c_gapsm,c_ordered,c_outputname)
