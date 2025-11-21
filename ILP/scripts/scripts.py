import random
import numpy as np










def gc_dfs(w,k,order,kmer,minorder,lev):
    n_kmers = 2**k
    n_kmers_mask = n_kmers-1
    return gc_dfs_prev_could_be_min(w,n_kmers_mask,order,kmer,minorder,lev)


# if here then preference is min so far
def gc_dfs_prev_could_be_min(w,n_kmers_mask,order,kmer,minorder,lev): #auxiliary#    dfs for prob_gc
    if lev == w: #leaf reached
        return 1
    else:
        nextkm = (kmer<<1) & n_kmers_mask # recursive call for children
        if order[kmer] >= minorder:
            return gc_dfs_prev_could_be_min(w,n_kmers_mask,order,nextkm,minorder=minorder,lev=lev+1) + gc_dfs_prev_could_be_min(w,n_kmers_mask,order,nextkm+1,minorder=minorder,lev=lev+1)
        else:    
            return gc_dfs_pref_is_not_min(w,n_kmers_mask,order,nextkm,minorder=order[kmer],lev=lev+1) + gc_dfs_pref_is_not_min(w,n_kmers_mask,order,nextkm+1,minorder=order[kmer],lev=lev+1)
        
        



def gc_dfs_pref_is_not_min(w,n_kmers_mask,order,kmer,minorder,lev):
    if lev == w: #leaf reached
        return order[kmer]<minorder # =1 if the leaf represents gamechanger
    else:
        nextkm = (kmer<<1) & n_kmers_mask # recursive call for children
        if order[kmer] >= minorder:
            return gc_dfs_pref_is_not_min(w,n_kmers_mask,order,nextkm,minorder=minorder,lev=lev+1) + gc_dfs_pref_is_not_min(w,n_kmers_mask,order,nextkm+1,minorder=minorder,lev=lev+1)
        else:    
            return gc_dfs_pref_is_not_min(w,n_kmers_mask,order,nextkm,minorder=order[kmer],lev=lev+1) + gc_dfs_pref_is_not_min(w,n_kmers_mask,order,nextkm+1,minorder=order[kmer],lev=lev+1)
        
def prob_gc(w,k,order): #fast count of gamechangers for an order of 2-ary k-mers
    total = 0
    n_kmers = 2**k
    n_kmers_mask = n_kmers-1
    for pref in range(n_kmers):
        nextkmer = (pref<<1) & n_kmers_mask
        total += gc_dfs(w,k,order,nextkmer,minorder=order[pref], lev=1) + gc_dfs(w,k,order,nextkmer+1, minorder=order[pref],lev=1)
    return total




def gc_dfs_function(w,k,function,kmer,minorder,lev):
    n_kmers = 2**k
    n_kmers_mask = n_kmers-1
    return gc_dfs_prev_could_be_min_function(w,k,n_kmers_mask,function,kmer,minorder,lev)


# if here then preference is min so far
def gc_dfs_prev_could_be_min_function(w,k,n_kmers_mask,function,kmer,minorder,lev): #auxiliary#    dfs for prob_gc
    if lev == w: #leaf reached
        return 1
    else:
        nextkm = (kmer<<1) & n_kmers_mask # recursive call for children
        order_kmer = function(kmer,w,k)
        if order_kmer >= minorder:
            return gc_dfs_prev_could_be_min_function(w,k,n_kmers_mask,function,nextkm,minorder=minorder,lev=lev+1) + gc_dfs_prev_could_be_min_function(w,k,n_kmers_mask,function,nextkm+1,minorder=minorder,lev=lev+1)
        else:    
            return gc_dfs_pref_is_not_min_function(w,k,n_kmers_mask,function,nextkm,minorder=order_kmer,lev=lev+1) + gc_dfs_pref_is_not_min_function(w,k,n_kmers_mask,function,nextkm+1,minorder=order_kmer,lev=lev+1)
        
        



def gc_dfs_pref_is_not_min_function(w,k,n_kmers_mask,function,kmer,minorder,lev):
    order_kmer = function(kmer,w,k)
    if lev == w: #leaf reached
        return order_kmer<minorder # =1 if the leaf represents gamechanger
    else:
        nextkm = (kmer<<1) & n_kmers_mask # recursive call for children
        if order_kmer >= minorder:
            return gc_dfs_pref_is_not_min_function(w,k,n_kmers_mask,function,nextkm,minorder=minorder,lev=lev+1) + gc_dfs_pref_is_not_min_function(w,k,n_kmers_mask,function,nextkm+1,minorder=minorder,lev=lev+1)
        else:    
            return gc_dfs_pref_is_not_min_function(w,k,n_kmers_mask,function,nextkm,minorder=order_kmer,lev=lev+1) + gc_dfs_pref_is_not_min_function(w,k,n_kmers_mask,function,nextkm+1,minorder=order_kmer,lev=lev+1)
        

def prob_gc_function(w,k,function): #fast count of gamechangers for an order of 2-ary k-mers
    total = 0
    n_kmers = 2**k
    n_kmers_mask = n_kmers-1
    for pref in range(n_kmers):
        nextkmer = (pref<<1) & n_kmers_mask
        order_pref = function(pref,w,k)
        total += gc_dfs_function(w,k,function,nextkmer,minorder=order_pref, lev=1) + gc_dfs_function(w,k,function,nextkmer+1, minorder=order_pref,lev=1)
    return total

def random_density_function(w, k, function):
    # create 10,000 w+k long binary strings
    strings = np.random.randint(0, 2, (10000, w + k))
    gc_count = 0
    for i in range(10000):
        string = int(''.join(map(str, strings[i])), 2)
        first_kmer = string & (2**k - 1)
        rank_first = function(first_kmer, w, k)
        last_kmer = string >> w
        rank_last = function(last_kmer, w, k)
        can_be_prefix = True
        can_be_suffix = True

        for j in range(1, w):
            kmer = (string >> (w - j)) & (2**k - 1)
            rank = function(kmer, w, k)
            if rank < rank_first:
                can_be_prefix = False
            if rank <= rank_last:
                can_be_suffix = False
            if not can_be_prefix and not can_be_suffix:
                break

        if can_be_prefix or can_be_suffix:
            gc_count += 1

    return gc_count / 10000
        



def rightdfs_wrapper(kmer,kmers,lev,w,k,gc, nongc):
    n_kmers = 2**k
    n_kmers_mask = n_kmers-1

    if lev == 0:
        return rightdfs_lev_zero(kmer,kmers,w,n_kmers_mask,gc, nongc)
    else:
        lev_left = w - lev
        return rightdfs(kmer,kmers,lev_left,n_kmers_mask,gc, nongc)            



def rightdfs(kmer,kmers,lev,w,n_kmers_mask,gc, nongc): #auxiliary#   recursively process a block of occurrences of kmer
    if nongc[kmer] == 0 and gc[kmer] == 0:
        return 0 #  stop as no more live windows contain kmer
    else:
        flag = kmer not in kmers
        if lev == w:
            gc[kmer] -= flag
            return 1
        else:
            if flag:
                kmers.add(kmer)
            nextkmer = (2*kmer) & n_kmers_mask
            sub = rightdfs(nextkmer,kmers,lev+1,w,n_kmers_mask,gc, nongc) + rightdfs(nextkmer+1,kmers,lev+1,w,n_kmers_mask,gc, nongc)
            if flag:
                if lev>0:
                    nongc[kmer] -= sub
                else:
                    gc[kmer] -= sub
                kmers.remove(kmer)
            return sub


def rightdfs_lev_zero(kmer,kmers,w,n_kmers_mask,gc, nongc):
    if nongc[kmer] == 0 and gc[kmer] == 0:
        return 0 #  stop as no more live windows contain kmer
    else:
        flag = kmer not in kmers
        nextkmer = (2*kmer) & n_kmers_mask

        if flag:
            kmers.add(kmer)
            sub = rightdfs(nextkmer,kmers,w-1,n_kmers_mask,gc, nongc) + rightdfs(nextkmer+1,kmers,w-1,n_kmers_mask,gc, nongc)
            gc[kmer] -= sub
            kmers.remove(kmer)
        else:
            sub = rightdfs(nextkmer,kmers,w-1,n_kmers_mask,gc, nongc) + rightdfs(nextkmer+1,kmers,w-1,n_kmers_mask,gc, nongc)

        return sub
    
def rightdfs(kmer,kmers,lev_left,n_kmers_mask,gc, nongc): #When lev left is 0 it means level = w
    if nongc[kmer] == 0 and gc[kmer] == 0:
        return 0 #  stop as no more live windows contain kmer
    else:
        flag = kmer not in kmers
        if lev_left == 0:
            gc[kmer] -= flag
            return 1
        else:
            nextkmer = (2*kmer) & n_kmers_mask
            if flag:
                kmers.add(kmer)
                sub = rightdfs(nextkmer,kmers,lev_left-1,n_kmers_mask,gc, nongc) + rightdfs(nextkmer+1,kmers,lev_left-1,n_kmers_mask,gc, nongc)
                nongc[kmer] -= sub
                kmers.remove(kmer)
            else:
                sub = rightdfs(nextkmer,kmers,lev_left-1,n_kmers_mask,gc, nongc) + rightdfs(nextkmer+1,kmers,lev_left-1,n_kmers_mask,gc, nongc)

            return sub

