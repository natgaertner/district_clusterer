import csv, difflib,operator,re
from pylab import hist, savefig,plot
import numpy as np
from itertools import combinations
from collections import defaultdict

class Cluster():
    def __init__(self, members, max_distances, min_distances):
        self.members = members
        self.max_distances = max_distances
        self.min_distances = min_distances

    def merge(self,cluster):
        max_distances = dict(map(lambda (sk,sv): max(((sk,sv),(sk,cluster.max_distances[sk])), key=lambda kv: kv[1]),((k,v) for k,v in self.max_distances.iteritems() if k != cluster)))
        min_distances = dict(map(lambda (sk,sv): max(((sk,sv),(sk,cluster.min_distances[sk])), key=lambda kv: kv[1]),((k,v) for k,v in self.min_distances.iteritems() if k != cluster)))
        return Cluster(self.members.extend(cluster.members),max_distances,min_distances)

def two_d_slice(array,imin,imax,jmin,jmax):
    return map(lambda inner: inner[jmin,jmax],array[imin:imax])

def two_d_min(array):
    outer_min,outer_mini = imin(array,key=lambda k:min(k))
    inner_min,inner_mini = imin(outer_min)
    return inner_min,outer_mini,inner_mini

def imin(array, key=lambda k:k):
    m = float('inf')
    mi = -1
    for i in xrange(len(array)):
        if key(array[i]) < m:
            m = key(array[i])
            mi = i
    return array[mi],mi

def two_d_max(array):
    outer_max,outer_maxi = imax(array,key=lambda k:max(k))
    inner_max,inner_maxi = imax(outer_max)
    return inner_max,outer_maxi,inner_maxi

def imax(array, key=lambda k:k):
    m = float('-inf')
    mi = -1
    for i in xrange(len(array)):
        if key(array[i]) > m:
            m = key(array[i])
            mi = i
    return array[mi],mi

def cluster_distance(cluster_one, cluster_two, ranks_one, ranks_two):
    return min(map(lambda (c1,r1):min(map(lambda (c2,r2): distance_metric(c1,c2,r1,r2),zip(cluster_two,ranks_two))),zip(cluster_one,ranks_one)))

def distance_metric(c1,c2,r1,r2):
    s = difflib.SequenceMatcher(None, c1,c2)
    match_blocks = s.get_matching_blocks()
    return sum(map(lambda b: sum(r1[b.a:b.a+b.size]) + sum(r2[b.b:b.b+b.size]),match_blocks))/float(sum(r1)+sum(r2))

def merge_clusters(cluster_one,cluster_two,ranks_one,ranks_two, min_block_size=1):
    ret_ranks_one = list(ranks_one)
    ret_ranks_two = list(ranks_two)
    for i in xrange(len(cluster_one)):
        c1 = cluster_one[i]
        r1 = ret_ranks_one[i]
        for j in xrange(len(cluster_two)):
            c2 = cluster_two[j]
            r2 = ret_ranks_two[j]
            r1,r2=merge_member(c1,c2,r1,r2, min_block_size)
            ret_ranks_two[j] = r2
        ret_ranks_one[i] = r1
    return cluster_one + cluster_two,tuple(ret_ranks_one + ret_ranks_two)

def merge_member(c1,c2,r1,r2, min_block_size=1):
    s = difflib.SequenceMatcher(None,c1,c2)
    matching_blocks = s.get_matching_blocks()
    retr1 = list(r1)
    retr2 = list(r2)
    for mb in matching_blocks:
        if mb.size >= min_block_size:
            for i in xrange(mb.a,mb.a+mb.size):
                retr1[i]+=1
            for i in xrange(mb.b,mb.b+mb.size):
                retr2[i]+=1
    return retr1,retr2

def get_fit_data(data):
    #X = np.arange(data.size)
    sorted_data = np.sort(data)
    u = float(sum(data))/len(data)
    s = float(sum(map(lambda x: (x-u)**2,data)))/len(data)
    print s
    if s == 0:
        def pdf(t):
            try:
                return map(lambda t_elem: 1 if t_elem == u else 0,t)
            except TypeError:
                return 1 if t == u else 0
        return pdf
    return lambda t: np.exp(-(t-u)**2/(2*s))/np.sqrt(2*np.pi*s)

def get_fit(u,s):
    if s == 0:
        def pdf(t):
            try:
                return map(lambda t_elem: 1 if t_elem == u else 0,t)
            except TypeError:
                return 1 if t == u else 0
        return pdf
    return lambda t: np.exp(-(t-u)**2/(2*s))/np.sqrt(2*np.pi*s)

ranks_flat = []
def calc_initial_ranks(clusters):
    merged = clusters[0]
    ranks = ([1]*len(merged[0]),)
    for c in clusters[1:]:
        merged, ranks = merge_clusters(merged,c,ranks,([1]*len(c[0]),), 3)
    global ranks_flat
    ranks_flat.extend([i for r in ranks for i in r])
    pdf = get_fit(125,15000)
    pdf_normed = lambda x:pdf(x)/pdf(125)
    return map(lambda r:(map(lambda i:pdf_normed(i),r),), ranks)

def plot_ranks():
    global ranks_flat
    hist(ranks_flat)
    savefig('hist.png')

class node(tuple):
    def __init__(self,*args,**kwargs):
        super(tuple,self).__init__(*args,**kwargs)
        self.children = set()
        self.count = 0

    def __hash__(self):
        return super(tuple,self).__hash__()

class passdict(dict):
    def __missing__(self,key):
        n = node(key)
        self[key] = n
        return n

words_patt = r'\b(?:\w|\s)+\b'
def is_child(parent, child):
    p_string = ' '.join(parent)
    p_string = re.sub(r'\s*\\b\(\?\:\\w\|\\s\)\+\\b\s*',words_patt,p_string)
    p_string = p_string.replace('\x08','\\b')
    pattern = re.compile(r'^'+p_string+r'$')
    return pattern.match(' '.join(filter(lambda t: t!=words_patt,child)))
    #return ' '.join(filter(lambda t: t != words_patt,parent)) in ' '.join(filter(lambda t: t != words_patt,child))

def classify_patterns(districts):
    #pattern_dict = defaultdict(lambda:0,{})
    pattern_dict = passdict()
    top_level = set()
    district_patterns = defaultdict(lambda:defaultdict(lambda:[],{}),{})
    district_pattern_counts = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:0,{}),{}),{})
    for d in districts:
        d_list = re.split(r'\s+',d)
        comb_dict = {}
        for r in range(len(d_list)):
            comb = combinations(range(len(d_list)),r+1)
            for c in comb:
                just_wrote=False
                if c[0] > 0:
                    patt = [words_patt]
                else:
                    patt = []
                for i in range(len(d_list)):
                    if i not in c and just_wrote:
                        patt.append(words_patt)
                        just_wrote=False
                    elif i in c and just_wrote:
                        patt.append(d_list[i])
                    elif i in c:
                        patt.append(d_list[i])
                        just_wrote=True
                patt = tuple(patt)
                comb_dict[c] = patt
                pattern_dict[patt].count+=1
                if r > 0:
                    for i in c:
                        c_prime = list(c)
                        c_prime.remove(i)
                        patt_prime = comb_dict[tuple(c_prime)]
                        pattern_dict[patt_prime].children.add(pattern_dict[patt])
                else:
                    top_level.add(pattern_dict[patt])
                district_patterns[d][r+1].append(patt)
    for d,dr in district_patterns.iteritems():
        for r,patts in dr.iteritems():
            for patt in patts:
                district_pattern_counts[d][r][patt] = pattern_dict[patt]
    test = dict((k,v) for k,v in pattern_dict.iteritems() if not is_child(pattern_dict[('MA', 'State', 'House', words_patt, 'District')],v) and is_child(pattern_dict[('MA','State',words_patt,'District')],v))
    test2 = map(lambda t: (t,t.count),test.values())
    test2.sort(key=lambda t:t[1])
    import pdb;pdb.set_trace()
    return pattern_dict,district_patterns


ordinal_patt = re.compile(r'\b(?P<number>\d+)(?:[sS][tT]|[rR][dD]|[nN][dD]|[tT][hH])\b')
def clean_ordinal(string):
    string_list = re.split(r'\s+',string)
    for i in xrange(len(string_list)):
        m = ordinal_patt.match(string_list[i])
        if m:
            string_list[i] = m.groupdict()['number']
    return ' '.join(string_list)

def clean_punc(string):
    return re.sub(r'[^A-Za-z0-9_ ]','',string)

def remove_nums(string):
    return re.sub(r'\d','',string)

def cluster_districts(districts, similarity_limit=.9, fix_ordinals=True, clean_punctuation=True, match_without_nums=False):
    districts = tuple(set(districts))
    if fix_ordinals:
        districts = tuple(map(lambda d: clean_ordinal(d),districts))
    if clean_punctuation:
        districts = tuple(map(lambda d: clean_punc(d),districts))
    clusters_track = map(lambda d: (d,),districts)
    if match_without_nums:
        districts = tuple(map(lambda d: remove_nums(d),districts))
        clusters = map(lambda d: (d,),districts)
    else:
        clusters = clusters_track
        #prox_array = map(lambda d: map(lambda innerd: difflib.SequenceMatcher(None,d,innerd).ratio(),districts_no_nums),districts_no_nums)
    #ranks = map(lambda c: ([1]*len(c[0]),),clusters)
    classify_patterns(districts)
    ranks = calc_initial_ranks(clusters)
    #prox_array = map(lambda d: map(lambda innerd: difflib.SequenceMatcher(None,d,innerd).ratio(),districts),districts)
    prox_array = map(lambda (d,r): map(lambda (innerd,innerr): cluster_distance(d,innerd,r,innerr),zip(clusters,ranks)),zip(clusters,ranks))
    for i in xrange(len(prox_array)):
        prox_array[i][i]=-1
    #import pdb;pdb.set_trace()
    m,i,j = two_d_max(prox_array)
    while m > similarity_limit:
        #print len(prox_array)
        #print len(clusters)
        #print m,i,j
        mi = min(i,j)
        mj = max(i,j)
        clusterj = clusters.pop(mj)
        clusteri = clusters.pop(mi)
        #rankj = ranks.pop(mj)
        #ranki = ranks.pop(mi)
        new_track_cluster = clusters_track.pop(mj) + clusters_track.pop(mi)
        #new_cluster,new_rank = merge_clusters(clusteri,clusterj,ranki,rankj)
        new_cluster = clusteri+clusterj
        proxj = prox_array.pop(mj)
        proxi = prox_array.pop(mi)
        prox_array = map(lambda p: p[0:mi]+p[mi+1:mj]+p[mj+1:],prox_array)
        proxi.pop(mj)
        proxi.pop(mi)
        proxj.pop(mj)
        proxj.pop(mi)
        prox = map(lambda (pi,pj):min(pi,pj),zip(proxi,proxj))
        #prox = map(lambda (c,r):cluster_distance(c,new_cluster,r,new_rank),zip(clusters,ranks))
        map(lambda i: prox_array[i].append(prox[i]),xrange(len(prox)))
        prox.append(-1)
        prox_array.append(prox)
        clusters.append(new_cluster)
        clusters_track.append(new_track_cluster)
        #ranks.append(new_rank)
        m,i,j = two_d_max(prox_array)
    return clusters_track

def detect_cluster_enumeration(cluster,min_block_size=1):
    enumerations = []
    match_strings = []
    for c in cluster:
        match_blocks = map(lambda oc: difflib.SequenceMatcher(None, ' '+c,' '+oc).get_matching_blocks(),cluster)
        #print map(lambda mb: map(lambda j: c[mb[j].a+mb[j].size:mb[j+1].a+1],xrange(len(mb)-1)),match_blocks)
        #import pdb;pdb.set_trace()
        poss_enum = max(map(lambda mb: reduce(operator.add,map(lambda j: (' '+c)[mb[j].a+mb[j].size:mb[j+1].a+1],xrange(len(mb)-1))),match_blocks),key=lambda m:len(m))
        enumerations.append(poss_enum)
        match_strings.append(min(map(lambda mb: reduce(operator.add,map(lambda b:(' ' +c)[b.a:b.a+b.size],filter(lambda b:b.size >= min_block_size,mb))),match_blocks),key=lambda m:len(m)))
        match_string = min(match_strings,key=lambda ms: len(ms))
    return enumerations,match_string

def guess_district(dist):
    dist = dist.lower()
    cong = ['congress','u.s. rep','us rep', 'u. s. rep', 'u s rep', 'united states rep']
    state_senate = ['state senat']
    state_rep = ['state house','state rep','legislat', 'state assembly']
    president = ['president']
    senate = ['us senat','u s senat','u.s. senat','u. s. senat','united states senat']
    if any(map(lambda d: d in dist,cong)):
        return 'congressional_district'
    if any(map(lambda d: d in dist,state_senate)):
        return 'state_senate_district'
    if any(map(lambda d: d in dist,state_rep)):
        return 'state_rep_district'
    if any(map(lambda d: d in dist,president)) or any(map(lambda d: d in dist,senate)):
        return 'state'
    else:
        return None

if __name__=='__main__':
    with open('dist_names.py','w') as dist_names:
        dist_names.write('{\n')
        #states = ['Arizona','Colorado','Georgia','New Mexico','Texas','California','Florida','Illinois','Pennsylvania','Massachusetts']
        states=['Massachusetts']
        for state in states:
            dist_names.write('{state}:{{'.format(state=repr(state)))
            with open('districts/{state}_districts.csv'.format(state=state),'w') as district_file,open('clusters/{state}_clusters'.format(state=state),'w') as cluster_file:
                offices = map(lambda r: r['office'],csv.DictReader(open('csvs/{state}.csv'.format(state=state))))
                clusters = cluster_districts(offices,clean_punctuation=True,match_without_nums=True,similarity_limit=.8)
                d= {}
                for c in clusters:
                    cluster_file.write(str(c))
                    cluster_file.write('\n')
                    enums, ms = detect_cluster_enumeration(c,3)
                    d.update(dict(zip(c,zip([ms]*len(enums),enums))))
                    dist_names.write("{ms}:".format(ms=repr(ms)))
                    dist_guess = guess_district(ms)
                    dist_names.write(repr(dist_guess)+",")
                    dist_names.write('\n')
                    cluster_file.write(str(enums))
                    cluster_file.write('\n')
                    cluster_file.write(ms)
                    cluster_file.write('\n\n')
                district_file.write(repr(d)+'\n')
                dist_names.write('},\n')
        dist_names.write('}')
        plot_ranks()
