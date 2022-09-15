import pandas as pd
import numpy as np
import random

import warnings
warnings.filterwarnings('ignore')

import math
import pickle

import pysam 
import pyBigWig

import itertools

from Bio.Seq import reverse_complement
from Bio.SeqIO.FastaIO import SimpleFastaParser

import joblib
from joblib import Parallel, delayed

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from denseweight import DenseWeight


def Fasta2Seqs(fileloc, idents = False, first = None): 
    
    seqs = []
    ids = []
    with open(fileloc) as fasta_file:
        for identifier, seq in SimpleFastaParser(fasta_file):
            seqs.append(seq.upper())
            ids.append(identifier)
            if first is not None and len(seqs) >= first: break
        
    return seqs, ids if idents == True else seqs

def CurrentSelect(Current, BS_ids, select_BS_ids): 
    idx = [BS_ids.index(z) for z in select_BS_ids if z in BS_ids] 
    return [Current[i] for i in idx]

def ListLengths(x): 
    return [len(y) for y in x] 

def Rounder(x, base=5):
    return base * np.round_(x/base).astype(int)

def Bam2Current(fileloc, BS_sizes, BS_ids, select_BS_ids = None, newreso = 20, resomode = np.sum, stranded = True, 
              read_length = None, paired = False, make_index = False, dtype = None): 
    
    #From 22.08.29.Bam2Current
    #When paired is True, the mates are mapped only by the first mate. Reads with unmapped mates are mapped as single end reads. 
    
    #22.09.14. added dtype for memory. could be np.int or np.float with 8/16/32. 
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids 

    if make_index == True: pysam.index(fileloc)
    
    samfile = pysam.AlignmentFile(fileloc, 'rb')
    
    T = []
    for z in select_BS_ids: 
        ind = BS_ids.index(z)
        ls = BS_sizes[ind]
        m = np.zeros((ls, 2), dtype = int)
        for read in samfile.fetch(z): 
            rp, rl, rt = read.pos, read.rlen, read.tlen
            
            rx = read_length if read_length is not None else rl 
            
            if paired == True and rt != 0: 
                if '{0:012b}'.format(read.flag)[::-1][6] == '1': rx = abs(read.tlen)
                else: continue 
                
            rp, st, e1, e2 = (rp + rl + 1, 1, rx, 0) if '{0:012b}'.format(read.flag)[::-1][4] == '1' else (rp, 0, 0, rx)
            j,k = rp-e1, rp+e2
            if j < 0: j = 0
            m[j:k, st] += 1
        
        if stranded == False: m = np.sum(m, axis = 1).reshape(-1,1)
    
        if newreso is not None:
            newlength = newreso * (BS_sizes[ind] // newreso)
            m = resomode(m[:newlength].reshape(newlength//newreso, newreso, -1), axis = 1)
            
        if dtype is not None: m = m.astype(dtype) 
        
        T.append(m) 
        
    return T

def Bw2Current(fileloc, BS_sizes, BS_ids, select_BS_ids = None, newreso = 20, resomode = np.sum, window = None, dtype = None):
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids 
    
    bw =  pyBigWig.open(fileloc)
    
    T = []
    for z in select_BS_ids: 
        ind = BS_ids.index(z)
        m = np.repeat(0, BS_sizes[ind])
        if window is not None: 
            for x in bw.intervals(z):
                if x[1]-x[0] < window: m[x[0]:x[1]] = x[2]
        else: 
            for x in bw.intervals(z): m[x[0]:x[1]] = x[2]
        if newreso > 1: 
            newlength = newreso * (BS_sizes[ind] // newreso)
            m = resomode(m[:newlength].reshape(-1, newreso), axis=1)
        
        if dtype is not None: m = m.astype(dtype) 
        T.append(m.reshape(-1,1))
    
    return T

def Bg2Current(fileloc, BS_sizes, BS_ids, select_BS_ids = None, newreso = 20, resomode = np.sum, window = None, sep = '\t', skiprows=None, dtype = None):
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids
    
    bg = pd.read_csv(fileloc, sep=sep, header = None, skiprows=skiprows)
    
    T = []
    for z in select_BS_ids: 
        ind = BS_ids.index(z)
        m = np.repeat(0, BS_sizes[ind])
        v = bg[bg[0] == z] 
        if window is not None: 
            for iro, ro in v.iterrows(): 
                if ro[2]-ro[1] < window: m[ro[1]:ro[2]] = ro[3]
        else: 
            for iro, ro in v.iterrows(): m[ro[1]:ro[2]] = ro[3]
        if newreso > 1: 
            newlength = newreso * (BS_sizes[ind] // newreso)
            m = resomode(m[:newlength].reshape(-1, newreso), axis=1)
            
        if dtype is not None: m = m.astype(dtype) 
        T.append(m.reshape(-1,1))
    
    return T

def Currents2Current(SigCurrents): 
    return [np.hstack([SigCurrents[T][z] for T in range(len(SigCurrents))]) for z in range(len(SigCurrents[0]))]

def CurrentMerger(SigCurrent, mode = np.mean, dtype = None):
    
    ms = []
    for z in SigCurrent: 
        m = mode(z, axis = 1).reshape(-1,1)
        if dtype is not None: m = m.astype(dtype)
        ms.append(m) 
    
    return ms

def CurrentModifier(SigCurrent, modifier, mode = np.multiply, dtype = None): 
    #modifier is a list, should match len of SigCurrent. Could be a single value ex. [-1] 
    
    mos = []
    for z,m in zip(SigCurrent, modifier):
        mo = mode(z, m)
        if dtype is not None: mo = mo.astype(dtype)
        mos.append(mo)
    
    return mos

def CurrentTransformer(SigCurrent, minmax = (0, None), standardize = False, dtype = None): 
    
    AT = [] 
    for t in range(SigCurrent[0].shape[1]): 
        T = [z[:,t] for z in SigCurrent]
        allvals = np.concatenate(T, axis = 0).reshape(-1)

        if standardize == True: 
            u = np.mean(allvals)
            w = np.std(allvals) 
            AT.append([((z-u)/w).astype(dtype).reshape(-1,1) for z in T]) 
        
        else: 
            mi, ma = minmax
            if mi is None: mi = np.min(allvals)
            if ma is None: ma = np.max(allvals) 
            mami = ma - mi
            AT.append([((z-mi)/mami).astype(dtype).reshape(-1,1) for z in T]) 
        
        
    return Currents2Current(AT)

def CurrentDtype(Current, dtype): 
    return [z.astype(dtype) for z in Current]

def CurrentInterpolator(SigCurrent, reso, window = None, threshold = 0, threshsign = np.greater): 
    
    AT = [] 
    for t in range(SigCurrent[0].shape[1]): 
        T = [z[:,t] for z in SigCurrent]
        nT = []
        for z in T: 
            lz = len(z)
            zpt = np.where(z > threshold)[0]
            zpt_vals = z[zpt]
            zI = np.interp(range(lz), zpt, zpt_vals)
            m = np.repeat(1, lz)

            if window is not None: 
                window = window//reso
                zpt = np.concatenate([[0], zpt , [lz]])
                for o in range(len(zpt)-1): 
                    oA,oB = zpt[o], zpt[o+1]
                    if oB-oA > window: 
                        m[oA:oB] = z[oA:oB]
            
            nT.append((zI*m).reshape(-1,1))
            
        if dtype is not None: nt = CurrentDtype(nT, dtype) 
        AT.append(nT) 
        
    return Currents2Current(AT)  

def LowerResCurrent(SigCurrent, reso, newreso, resomode = np.mean, dtype = None): 
    
    resoratio = newreso // reso
    AT = [] 
    for t in range(SigCurrent[0].shape[1]): 
        T = [z[:,t] for z in SigCurrent]
        nT = []
        for z in T: 
            newlength = resoratio * (len(z)//resoratio)
            nT.append(resomode(z[:newlength].reshape(-1, resoratio), axis = 1).reshape(-1,1))
        
        if dtype is not None: nt = CurrentDtype(nT, dtype) 
        AT.append(nT)
    
    return Currents2Current(AT)

def CurrentExpander(SigCurrent, reso, BS_size, BS_ids, select_BS_ids = None, filler = 0): 
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids
    
    nT = []
    for iz, z in enumerate(select_BS_ids):
        ind = BS_ids.index(z)
        lg = BS_size[ind]
        repd = np.repeat(SigCurrent[iz].reshape(-1), reso)
        dy = SigCurrent[iz].dtype.type
        nT.append(np.append(repd, np.repeat(filler,lg-len(repd)).astype(dy)).reshape(-1,1))
        
    return nT

def BinCurrentInverter(BinCurrent): 
    return [1-b for b in BinCurrent] 

def CurrentWindower(SigCurrent, reso = 1, mode = np.mean, window = None, center = None, extend = None, dtype = None): 
    
    window = 0 if window is None else window//reso 
    window = 1 if window == 0 else window
    
    TM = []
    for z in SigCurrent: 
        lz = len(z) 
        r = np.lib.stride_tricks.sliding_window_view(z, (window, z.shape[1])) 
        t = mode(mode(r, -1),-1).reshape(-1)
        if dtype is not None: t = t.astype(dtype)
        dy = t.dtype.type
        
        if center is not None: t = np.concatenate([np.repeat(center, window // 2).astype(dy), t])
        if extend is not None: t = np.concatenate([t, np.repeat(extend, lz - len(t)).astype(dy)])
        
        TM.append(t.reshape(-1,1))
    
    return TM 

def CurrentThresholder(SigCurrent, threshold = 0, threshsign = np.greater, binary = False): 
    
    TT = [(threshsign(z, threshold)*1).astype(np.int8) for z in SigCurrent]
    if binary is not True: TT = [z*t for z,t in zip(SigCurrent, TT)]
    
    return TT

def CurrentPadder(SigCurrent, reso = 1, pad = 10, shift = 0, value = 1, base = 0, dtype = None): 
    
    pad, shift = (0 if x is None else x//reso for x in [pad, shift])
    
    TB = []
    for z in SigCurrent: 
        lz = len(z)
        m = np.repeat(base, lz) if dtype is None else np.repeat(base, lz).astype(dtype)
        v = np.where(z == value)[0]+shift
        pdx = np.unique(np.concatenate([np.arange(w - pad, w + pad) for w in (np.where(z == value)[0]+shift)]))
        pdx = pdx[pdx<lz]
        m[pdx] = value 
        
        TB.append(m.reshape(-1,1))
    
    return TB

def CurrentFlagger(SigCurrent, reso = 5, window = 1000, stride = 50, flagmode = np.argmax, double = False): 
    #flagmode can be either np.argmax, np.argmin or any index searcher 
    
    window, stride = (0 if x is None else x//reso for x in [window,stride])
    window, stride = (1 if x == 0 else x for x in [window,stride])
    
    j = SigCurrent[0].shape[1]
    ms = []
    for z in SigCurrent: 
        lz = len(z)
        r = np.lib.stride_tricks.sliding_window_view(z, (window, j))
        
        idx = np.unravel_index([flagmode(x) for x in r[::stride]], (window,j))[0] + np.arange(0, lz-window+1, stride)
        #idx = np.unravel_index(flagmode(r.reshape(r.shape[0], -1), axis = 1), (window, j))[0] + np.arange(lz - window + 1)
        idx = np.unique(idx)
        
        m = np.repeat(0, lz).astype(np.int8).reshape(-1,1)
        m[idx] = 1
        
        if double == True: 
            
            y = z.copy()
            window2 = window//2
            
            for ix in idx: 
                
                g1, g2 = (np.zeros((window2, f)).astype(np.int8) for f in [j, 1])
                ix2 = np.unravel_index(flagmode(y[ix:ix+window2], axis = 1), (window2, j))[0]
                g1[ix2], g2[ix2] = 1, 1
                
                if ix+window2 > lz: 
                    g1, g2 = (x[:lz-ix] for x in [g1, g2])
            
                y[ix:ix+window2] = y[ix:ix+window2] * g1
                m[ix:ix+window2] = m[ix:ix+window2] * g2

        ms.append(m)
    
    return [m*z for m,z in zip(ms, SigCurrent)]


def Seqs2kmer(seqs, k = 6): 
    return [[seq[x:x+k] for x in np.arange(len(seq) - k + 1)] for seq in seqs]

def kmerFinder(seqs, kmer, RC = True, center = True, extend = True):
    lk = len(kmer)
    KF = [np.array([1 if x == kmer else 0 for x in Seqs2kmer([seq], k = lk)[0]]) for seq in seqs]
    
    if RC == True: 
        rKF = [np.array([1 if x == kmer else 0 for x in Seqs2kmer([reverse_complement(seq)], k = lk)[0]]) for seq in seqs]
        rKF = [x[::-1] for x in rKF] 
        KF = [np.max([u,w], axis = 0) for u,w in zip(KF, rKF)]
    
    if center == True: KF = [np.concatenate([np.repeat(0, lk // 2), x]) for x in KF]
    if extend == True: KF = [np.concatenate([x, np.repeat(0, len(seqs[ix]) - len(x))]) for ix,x in enumerate(KF)]
    
    return [x.reshape(-1,1) for x in KF]

def Seqs2OHE(seqs, expand_dim = None, customdict = None):
    
    OHEdict = {'A': [1,0,0,0],'C': [0,1,0,0],'G': [0,0,1,0],'T': [0,0,0,1],
               'N': [0.25,0.25,0.25,0.25],'R': [0.5,0,0.5,0], 'Y': [0,0.5,0,0.5], 
               'S': [0,0.5,0.5,0], 'W': [0.5,0,0,0.5], 'K': [0,0,0.5,0.5], 
               'M': [0.5,0.5,0,0], 'B': [0,0.33,0.33,0.33], 'D': [0.33,0,0.33,0.33], 
               'H': [0.33,0.33,0,0.33], 'V': [0.33,0.33,0.33,0] }
    if customdict is not None: OHEdict = customdict
    
    OHEtotal = []
    for seq in seqs:
        seq = seq.upper()
        OHEseq = np.array([OHEdict[nuc] for nuc in seq]).astype(np.float16)
        OHEtotal.append(OHEseq)
    OHEtotal = np.squeeze(np.stack(OHEtotal))
    if expand_dim is not None: OHEtotal = np.expand_dims(OHEtotal, axis = expand_dim)
    
    return OHEtotal

def Seqs2OHERC(seqs, expand_dim = None, customdict = None):

    seqs_rc = [reverse_complement(x) for x in seqs]
    
    OHEd, OHEd_rc = (Seqs2OHE(x, customdict = customdict) for x in [seqs, seqs_rc])
    OHE_comb = np.concatenate((OHEd, OHEd_rc), axis = -1)
    if expand_dim is not None: OHE_comb = np.expand_dims(OHE_comb, axis = expand_dim)
    
    return OHE_comb

def CM2PWM(fileloc, pseudo = True, bkg = [0.25,0.25,0.25,0.25], base = 2): 
    
    CM = np.loadtxt(fileloc)
    
    if pseudo == True:
        CM = CM + 1 
    poz = []
    for p in CM: 
        tot = sum(p) 
        pos = []
        for i,n in enumerate(p): 
            pos.append(math.log((n/tot)/bkg[i], base))
        poz.append(pos) 
    
    return np.array(poz) 

def Seq2PWMScore(seq, PWM, pieces = None, RC = True, center = None, extend = None, dtype = None): 
    # introduce pieces for long seq 
    
    ls, lp = len(seq), len(PWM)
    
    pieces = 1 if pieces == None else pieces
    parts = np.linspace(0, ls - lp + 1, pieces + 1, dtype = int)
    
    scores = []
    for i in range(pieces): 
        sp = seq[parts[i]:parts[i+1] + lp - 1]
        s = Seqs2OHE(sp)
        r = np.lib.stride_tricks.sliding_window_view(s, (lp, 4))
        s1 = np.tensordot(r, PWM, axes=((2,3),(0,1))).reshape(-1)
        del s,r
        if RC == True: 
            s_rc = Seqs2OHE(reverse_complement(sp))
            r_rc = np.lib.stride_tricks.sliding_window_view(s_rc, (lp, 4))
            s2 = np.tensordot(r_rc, PWM, axes=((2,3),(0,1))).reshape(-1)
            del s_rc,r_rc
            s1 = np.max([s1, s2[::-1]], axis = 0)
        
        if dtype is not None: s1 = s1.astype(dtype) 

        scores.append(s1)
    
    scores = np.concatenate(scores)
    dy = scores.dtype.type
    if center is not None: scores = np.concatenate([np.repeat(center, lp // 2).astype(dy), scores])
    if center is not None: scores = np.concatenate([scores, np.repeat(extend, ls - len(scores)).astype(dy)])
        
    return scores.reshape(-1,1)

def PWMScorer(PWM, BaseSeq, BS_ids, select_BS_ids = None, parrellel = None, pieces = None, RC = True, center = None, extend = None, dtype = None): 
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids
    
    lp = len(PWM)
    
    pieces = 1 if pieces == None else pieces
    
    PWMscores = []
    
    for t in select_BS_ids: 
        if t in BS_ids: 
            ind = BS_ids.index(t)
            g = BaseSeq[ind]
            ls = len(g)
            parts = np.linspace(0, ls - lp + 1, pieces + 1, dtype = int)
            if parrellel is not None: 
                delayed_funcs = [delayed(Seq2PWMScore)(g[parts[i]:parts[i+1] + lp - 1], PWM, pieces = 1, dtype = dtype) for i in range(pieces)]
                parallel_pool = Parallel(n_jobs=parrellel, backend = 'multiprocessing') #to make sure things are in order use multiprocessing
                scores = np.vstack(parallel_pool(delayed_funcs))
                dy = scores.dtype.type
                if center is not None: scores = np.vstack([np.repeat(center, lp // 2).astype(dy).reshape(-1,1), scores])
                if center is not None: scores = np.vstack([scores, np.repeat(extend, ls - len(scores)).astype(dy).reshape(-1,1)])
                PWMscores.append(scores) 
                
            else: PWMscores.append(Seq2PWMScore(g,PWM,pieces = pieces, RC = RC, extend = extend, center = center, dtype = dtype))
    
    return PWMscores

def Bed2Markers(fileloc, sep = '\t', header = None, skiprows = None, strandcol = None, addcols = None, addcols_names = None): 
    
    Markercols = ['subseq', 'loc', 'strand', 'size']
    if addcols is not None: 
        Markercols = Markercols + ['add' + str(x) for x in addcols] if addcols_names is None else Markercols + addcols_names
    px = pd.read_csv(fileloc, sep=sep, header = header,skiprows = skiprows)
    
    Markers = []
    for i in range(len(px)): 
        p = px.iloc[i]
        st = p[strandcol] if strandcol is not None else 0
        si = p[2] - p[1]
        k = (p[q] for q in addcols) if addcols is not None else []
        Markers.append([p[0], (p[1] + p[2]) // 2, st, si, *k])
    
    Markers = pd.DataFrame(Markers, columns = Markercols).sort_values(by=['subseq', 'loc', 'strand']).reset_index(drop = True) 
    
    return Markers

def MarkersFilter(Markers, SigCurrent, BS_ids, select_BS_ids = None, reso = 1, exact = False, Msizes = None, threshold = None, threshmode = np.sum, threshsign = np.greater): 
    
    #Msizes could be Markers[:, 3]
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids
    
    MarkersFilt = Markers[Markers.iloc[:, 0].isin(select_BS_ids)].reset_index(drop = True)
    
    si = np.repeat(1, len(Markers)) if Msizes is None else Msizes
    if len(Msizes) == 1: si = np.repeat(Msizes[0], len(Markers))
    
    pp = []
    for ir in MarkersFilt.index: 
        V = SigCurrent2Slice(SigCurrent, MarkersFilt.iloc[ir] , BS_ids, select_BS_ids, Vsize = si[ir], reso = reso, newreso = None, resomode = np.mean, exact = exact, Vpad = 2)
        if threshsign(threshmode(V, axis = 1), threshold) == True: pp.append(ir)
    
    MarkersFilt = MarkersFilt.iloc[pp].reset_index(drop = True) 
    
    return MarkersFilt

def Markers2BinCurrent(Markers, BS_size, BS_ids, select_BS_ids = None, Msizes = None, reso = 1, inverse = False, dtype = None): 
    #Msizes could be Markers[:, 3]
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids
    
    si = np.repeat(1, len(Markers)) if Msizes is None else Msizes
    if len(Msizes) == 1: si = np.repeat(Msizes[0], len(Markers))
    si = np.array([(s//reso)//2 for s in si])
    sp = np.array([1 if s == 0 else 0 for s in si])
    
    x,y = (0,1) if inverse == False else (1,0)
    
    ms = []
    for iz,z in enumerate(select_BS_ids): 
        ms.append(np.repeat(x, BS_size[BS_ids.index(z)]//reso))
        if dtype is not None: ms[iz] = ms[iz].astype(dtype)
    
    for ir, (z,r) in enumerate(zip(Markers.iloc[:, 0], Markers.iloc[:, 1])): 
        if z in select_BS_ids: 
            ind = select_BS_ids.index(z) 
            f = Rounder(r, reso) 
            k = si[ir] 
            ms[ind][f-k:f+k+sp[ir]] = y
    
    return [m.reshape(-1,1) for m in ms] 


def SigCurrent2Markers(SigCurrent, BS_ids, select_BS_ids = None, reso = 1, select = None, select_mode = None, mode_args = None): 
    
    if select is not None: 
        
        allvals = np.concatenate(SigCurrent, axis = 0).reshape(-1)
        nz_idx = np.where(allvals != 0)[0]
        nz_vals = allvals[nz_idx] 
        
        if select_mode is not None: 
            mode_idx = select_mode(nz_vals, select, **mode_args if mode_args is not None else None)
            sel_idx = nz_idx[mode_idx]
        
        else: sel_idx = np.random.choice(nz_idx, replace = False, size = select)
        
        pdx = IdxFlat2Current(SigCurrent, sel_idx)
    
    else: pdx = [np.where(z != 0) for z in SigCurrent]
    
    Markers = []
    for ip, p in enumerate(pdx): 
        lp = len(p)
        vals = (SigCurrent[ip][p]).reshape(-1)
        subseqs, strand, size = [np.repeat(x, lp) for x in [select_BS_ids[ip], 0, 1]]
        Markers.append(np.stack([subseqs, p*reso, strand, size, vals], axis = 1))
    Markers = np.vstack(Markers)
    
    d = pd.DataFrame(Markers, columns = ['subseq', 'loc', 'strand', 'size', 'val'])
    d[['loc', 'strand', 'size']] = d[['loc', 'strand', 'size']].astype(int)
    d[['val']] = d[['val']].astype(float)
    
    d = d.sort_values(by=['subseq', 'loc', 'strand']).reset_index(drop = True)
    
    return d    

def Top(vals, select, smallest = False): 
    
    k, m = (-1, np.greater_equal) if smallest == False else (1, np.less_equal)
    
    thresh = vals[np.argsort(vals)[::k]][select]
    
    idx = np.where(m(vals, thresh))[0]
    if len(idx) > select: idx = np.random.choice(idx, replace = False, size = select)
    
    return idx 

def Harpoon(vals, select, bins = 100):
    
    h1, h2 = np.histogram(vals, bins)
    div = [np.where((vals>=h2[i]) & (vals<=h2[i+1]))[0] for i in range(bins)]
    
    j = 0 
    for v in range(select): 
        for i in range(bins): 
            if h1[i] > 0: 
                j += 1
                h1[i] -= 1
        if j >= select: 
            w = v + 1
            break
    
    idx = np.hstack([d if len(d) < w else np.random.choice(d, replace = False, size = w) for d in div]) 
    if len(idx) > select: idx = np.random.choice(idx, replace = False, size = select)
    
    return idx

def RevDistro(vals, select, alpha = 1.0): 
    
    dw = DenseWeight(alpha=alpha)
    weights = dw.fit(vals)
    weights = weights / np.sum(weights)

    idx = np.random.choice(np.arange(len(vals)), size = select, replace = False, p = weights)
    
    return idx

def IdxFlat2Current(SigCurrent, flatidx): 
    
    lz = [0] + [len(z) for z in SigCurrent]
    rs = np.cumsum(lz)
    Currentidx = [flatidx[np.where((flatidx >= rs[i]) & (flatidx < rs[i+1]))[0]] - rs[i] for i in range(len(SigCurrent))]
    
    return Currentidx

def Txt2Markers(fileloc, addcols, sep = '\t', header = None, skiprows = None, addcols_names = None):
    
    Markercols = ['subseq', 'loc', 'strand', 'size']
    if addcols is not None: 
        Markercols = Markercols + ['add' + str(x) for x in addcols] if addcols_names is None else Markercols + addcols_names
    px = pd.read_csv(fileloc, sep=sep, header = header)
    
    Markers = []
    for i in range(len(px)): 
        p = px.iloc[i]
        k = (p[q] for q in addcols)
        Markers.append([p[0], 0, 0, 0, *k])
    
    Markers = pd.DataFrame(Markers, columns = Markercols).sort_values(by=['subseq']).reset_index(drop = True)
    
    return Markers

def AddTraits(Markers_A, Markers_B, Msizes_A = None, Msizes_B = None): 
    
    if Msizes_A is None and Msizes_B is None: 
        merged = pd.merge(Markers_A, Markers_B, how = 'left', on='subseq')
    
    else: 
    
        lR = len(Markers_A)
        si_A, si_B = (np.repeat(1, lR) if R is None else R for R in [Msizes_A, Msizes_B])
        si_A, si_B = (np.repeat(R[0], lR) if len(R) == 1 else R for R in [Msizes_A, Msizes_B])
        si_X = (si_A + si_B) // 2
        j,k = Markers_B['loc'] - si_X, Markers_B['loc'] + si_X
        
        attr = []
        for ir, r in enumerate(Markers_A): 
            at = Markers_B.loc[(Markers_B['subseq'] == r[0]) & (j <= r[1]) & (k >= r[1])]
            if len(at) > 0: attr.append(at.iloc[0].tolist())
            else: attr.append([])
        attr = pd.DataFrame(attr, columns = Markers_B.columns[4:])
        
        merged = pd.concat([Markers_A, attr])
        
    return merged 

def Markers4Packs(Markers, select_BS_ids = None, opposite = False, ends = None, BS_sizes = None, BS_ids = None): 
    
    if select_BS_ids is not None: Markers = Markers[Markers['subseq'].isin(select_BS_ids)].reset_index(drop = True)
    
    if ends is not None: 
        ss = Markers['subseq'].unique()
        goodidx = []
        for u in ss: 
            ind = BS_ids.index(u)
            le = BS_sizes[ind] - ends
            g = Markers[Markers['subseq'] == u]
            goodidx.append(g[(g['loc'] - g['size'] > ends) & (g['loc'] + g['size'] < le)].index)
        goodidx = np.concatenate(goodidx)
        Markers = Markers.iloc[goodidx]
        
    if opposite == True: 
        MarkersOpp = Markers.copy() 
        MarkersOpp['strand'] = 1 - MarkersOpp['strand']
        Markers = pd.concat([Markers, MarkersOpp]).sort_values(by=['subseq', 'loc', 'strand'])
        Markers = Markers.drop_duplicates().reset_index(drop = True)
    
    return Markers

def SeqCurrent2Slice(SeqCurrent, Marker, BS_ids, select_BS_ids, Vsize = 1000, stranded = False, seqmode = None): 
    #seqmode is either None or Seqs2OHE or Seqs2OHERC 
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids
    
    ind = select_BS_ids.index(Marker[0])
    ud = Vsize // 2
    r = Marker[1]
    V = SeqCurrent[ind][r-ud:r+ud]
    if stranded == True: 
        if Marker[2] > 0: V = reverse_complement(V)
    
    return seqmode([V]) if seqmode is not None else V 

def SeqCurrent2Pack(SeqCurrent, Markers, BS_ids, select_BS_ids = None, Vsize = 1000, seqmode = None, stranded = False): 
    #mode is either Seqs2OHE or Seqs2OHERC 
    #if RC is true, it will add the reverse complement of the pack at the end 
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids
    
    seqs = [] 
    for ir in Markers.index: 
        seqs.append(SeqCurrent2Slice(SeqCurrent, Markers.iloc[ir], BS_ids, select_BS_ids, Vsize = Vsize, stranded = stranded, seqmode = None))
    
    if seqmode is not None: seqs = seqmode(seqs) 
    
    return seqs

def SigCurrent2Slice(SigCurrent, Marker, BS_ids, select_BS_ids = None, Vsize = 1000, reso = 5, newreso = None, resomode = np.mean, 
                     exact = True, Vpad = 2, stranded = False, dtype = None): 
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids
    ind = select_BS_ids.index(Marker[0])

    ud = (Vsize//2) // reso
    p = 1 if ud == 0 else 0 
    r = Marker[1]
    j = Rounder(r,reso)
    c = j // reso
    
    if reso == 1: exact = False
    
    if exact == True:
        g = (Vpad*reso) + (r - j)
        V = SigCurrent[ind][c-ud-Vpad : c+ud+Vpad+p]
        V = np.stack([np.interp(np.arange(g, g+Vsize, reso), np.arange(0, ((ud+Vpad)*2*reso)+p, reso), b) for b in V.T], axis = 0)
    else: 
        V = SigCurrent[ind][j-ud:j+ud+p].T

    if stranded == True: 
        if Marker[2] > 0: V = np.flip(V, axis = -1)
    
    if newreso is not None: 
        resoratio = 1 if newreso is None else newreso // reso
        newlength = resoratio * ((V.shape[-1])//resoratio)
        V = resomode(V[:, :newlength].reshape(V.shape[0], -1, resoratio), axis = -1)
    
    if dtype is not None: V = V.astype(dtype)
    
    return V

def SigCurrent2Pack(SigCurrent, Markers, BS_ids, select_BS_ids = None, Vsize = 1000, reso = 5, newreso = None, resomode = np.mean, 
                    areas = None, exact = False, stranded = False, dtype = None): 
    #newreso needs to be divisiable by reso
    #areas is a list of areas to output instead of profile/single area
    
    select_BS_ids = BS_ids if select_BS_ids is None else select_BS_ids
    
    Pack = []
    for ir in Markers.index:
        Pack.append(SigCurrent2Slice(SigCurrent, Markers.iloc[ir] , BS_ids, select_BS_ids, Vsize = Vsize, reso = reso, newreso = None, 
                                     resomode = np.mean, exact = exact, Vpad = 2, stranded = stranded))

    Pack = np.stack(Pack, axis = 0)
    
    if newreso is not None: 
        resoratio = 1 if newreso is None else newreso // reso
        newlength = resoratio * ((Pack.shape[-1])//resoratio)
        Pack = resomode(Pack[:, :, :newlength].reshape(Pack.shape[0], Pack.shape[1], -1, resoratio), axis = -1)
    
    rr = reso if newreso is None else reso*resoratio
    if areas is not None: 
        ax = []
        c = Pack.shape[-1] // 2
        for a in areas: 
            q = (a//2)//rr
            ax.append(np.trapz(Pack[:, :, c-q:c+q], axis = -1))
        Pack = np.stack(ax, axis = -1)
    
    if dtype is not None: Pack = Pack.astype(dtype)
    
    return Pack

def LowerResProfilePack(ProfilePack, reso, newreso, resomode = np.mean, dtype = None): 
    
    resoratio = newreso // reso
    newlength = resoratio * ((ProfilePack.shape[-1])//resoratio)
    ProfilePack = mode(ProfilePack[:, :, :newlength].reshape(ProfilePack.shape[0], ProfilePack.shape[1], -1, resoratio), axis = -1)
    
    if dtype is not None: ProfilePack = ProfilePack.astype(dtype)
    
    return ProfilePack

def Traits2Pack(Markers, cols = None, vectorize = None, dtype = None): 
    
    x = Markers.iloc[:, 4:] if cols is None else Markers.iloc[:, cols]
    
    if vectorize is not None: 
        names = Markers.columns[vectorize] 
        x = pd.get_dummies(x, prefix = names)
    
    x = x.to_numpy()
    if dtype is not None: x = x.astype(dtype)
    
    return x

def PackShaper(Pack, expand_dim = None, squeeze = False, twoD = False, flatten = False): 
    
    if expand_dim is not None: Pack = np.expand_dims(Pack, axis = expand_dim)
    if squeeze == True: Pack = np.squeeze(Pack) 
    if twoD == True: Pack = np.reshape(Pack, (len(Pack), -1))
    if flatten == True: Pack = np.reshape(Pack, (-1)) 
    
    return Pack

def Packs2Pack(Packs): 
    return np.hstack(Packs)

def Splitter(X, numsplits = 3, proportions = [0.7, 0.3], random = False, cut = True): 
    
    #X is the len of Pack which is the same len of Markers 
    
    idx = np.arange(X)
    if random == True: np.random.shuffle(idx) 

    proportions = np.array(proportions) / np.sum(proportions)
    nums = np.round(proportions * X).astype(int)
    v = np.arange(len(proportions))
    
    SplitDict = {}
    for i in range(numsplits): 
        SplitDict[i] = {}
        np.random.shuffle(v)
        if cut == True: 
            c = np.random.choice(idx)
            idx = np.concatenate([idx[c:], idx[0:c]]) 
        counter = 0
        for e in v: 
            p = nums[e]
            SplitDict[i][e] = idx[counter:counter + p]
            counter += p
    
    Splits = [[SplitDict[x][i] for i in range(len(proportions))] for x in SplitDict.keys()]
    
    return Splits

def PackSplit(Pack, Split): 
    return [Pack[s] for s in Split]

def GridPlot(Arr, bounds = (None, None), 
                xlabels = None, labels = None, suptitle = None, 
                figsize = (15,5), fontsize = 10, cmap = 'viridis'): 
    
    fig, ax = plt.subplots(figsize = figsize)
    
    g = ax.imshow(Arr, vmin = bounds[0], vmax = bounds[1], cmap = cmap)
    
    if labels is not None: 
        ax.set_yticks(np.linspace(0, Arr.shape[0]-1, len(labels)))
        ax.set_yticklabels(labels, fontsize = fontsize) 
    if xlabels is not None: 
        ax.set_xticks(np.linspace(0, Arr.shape[1]-1, len(xlabels)))
        ax.set_xticklabels(xlabels, fontsize = fontsize)
        #plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    if suptitle is not None: ax.set_title(suptitle, size = fontsize, loc = 'left')
    
    return

def FilledLinePlot(Arr, bounds = (None, None), 
                xlabels = None, labels = None, suptitle = None, 
                figsize = (15,8), fontsize = 10, cmap = 'viridis'): 
    
    colors = plt.get_cmap(cmap, len(Arr))
    k = Arr.shape[1]
    
    sharey = True if bounds[0] is None and bounds[1] is None else False 
    
    fig, axs = plt.subplots(len(Arr), 1, figsize = figsize, constrained_layout=True, sharey = sharey)
    
    for i, d in enumerate(Arr): 
        p = axs if len(Arr) == 1 else axs[i]
        p.fill_between(np.arange(k), 0, d, color = colors(i))
        
        p.spines['right'].set_visible(False)
        p.spines['top'].set_visible(False)
        
        if sharey == False: p.set_ylim(bounds) 
        
        if xlabels is not None:  
            p.set_xticks(np.linspace(0, k-1, len(xlabels)), labels = [])
            if i == 0: p.set_xticklabels(xlabels, fontsize = fontsize)
        
        p.margins(x = 0)
        
        if labels is not None: p.set_title(labels[i], loc = 'left', fontsize = fontsize)
    
    plt.suptitle(suptitle, x = 0, size = fontsize)
    
    return 

def LinePlot(Arr, bounds = (None, None), 
                xlabels = None, labels = None, suptitle = None, 
                figsize = (15,8), fontsize = 10, cmap = 'viridis'): 
    
    colors = plt.get_cmap(cmap, len(Arr))
    k = Arr.shape[1]
    
    sharey = True if bounds[0] is None and bounds[1] is None else False 
    
    fig, axs = plt.subplots(len(Arr), 1, figsize = figsize, constrained_layout=True, sharey = sharey)
    
    for i, d in enumerate(Arr): 
        p = axs if len(Arr) == 1 else axs[i]
        p.plot(np.arange(k), d, color = colors(i))
        
        p.spines['right'].set_visible(False)
        p.spines['top'].set_visible(False)
        
        if sharey == False: p.set_ylim(bounds) 
        
        if xlabels is not None:  
            p.set_xticks(np.linspace(0, k-1, len(xlabels)), labels = [])
            if i == 0: p.set_xticklabels(xlabels, fontsize = fontsize)
        
        p.margins(x = 0)
        
        if labels is not None: p.set_title(labels[i], loc = 'left', fontsize = fontsize)
    
    plt.suptitle(suptitle, x = 0, size = fontsize)
    
    return 

def SigCurrent2Visual(SigCurrent, Marker, BS_ids, select_BS_ids = None, reso = 5, newreso = None, exact = False, 
                 Vsize = 2000, vismode = GridPlot, 
                 bounds = (None, None), xlabelsnum = 3, labels = None, suptitle = None, figsize = (15,5), fontsize = 10, cmap = 'viridis'):
    
    #Mode can be GridPlot or oor LinePlot or FilledLinePlot(for now) 
    
    V = SigCurrent2Slice(SigCurrent, Marker, BS_ids, select_BS_ids = select_BS_ids, Vsize = Vsize, reso = reso, newreso = newreso, exact = exact, resomode = np.mean, Vpad = 2)
    
    ud = Vsize // reso // 2
    x = Marker[1] if exact == True else Marker[1] // reso
    xlabels = np.linspace(x - (ud*reso), x + (ud*reso), xlabelsnum, dtype = int)
    
    vismode(V, bounds = bounds, xlabels = xlabels, labels = labels, suptitle = suptitle, figsize = figsize, fontsize = fontsize, cmap = cmap)

    return

def ProfilePack2Visual(ProfilePack, combomode = np.mean, 
                       Vsize = 1000, vismode = GridPlot, 
                       bounds = (None, None), xlabelsnum = 3, labels = None, suptitle = None, figsize = (15,5), fontsize = 10, cmap = 'viridis'):
    
    #ProfilePack format is (regions, Currents, profile) 
    #combomode can be np.mean or np.sum
    #Mode can be GridPlot or FilledLines(for now) 
    
    cV = combomode(ProfilePack, axis = 0) 
    xlabels = np.linspace(0, Vsize, xlabelsnum, dtype = int)
    vismode(cV, bounds = bounds, xlabels = xlabels, labels = labels, suptitle = suptitle, figsize = figsize, fontsize = fontsize, cmap = cmap)
