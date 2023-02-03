import re
import sys
import time

sys.path.insert(0,"/mnt/extraids/SDSC_NFS/linpei/bin/HiFiToolkit")
import hifiwrangler0

spot_to_hifi_sam = sys.argv[1]
spot_duplic_info = sys.argv[2]

Nmax = 0
if len(sys.argv) == 4:
    Nmax = int(sys.argv[3])

HiFi = {}
Spot_dic = {}
Spot_obj = {}
Spot_2_N = {}

# [1] https://www.geeksforgeeks.org/how-to-read-large-text-files-in-python/
# in [1] I found that this way of reading file (for example, lines 24 - 25) could be faster than readline() and input()
# I did not compare this strategy with Pandas.

with open(spot_to_hifi_sam) as file:
    for line in file:
        a = line.split("\t")
        for i in a:
            if i.find("AS:i:") != -1:
                # this field is the AS score. We will create an instance of class Read using the a[0]. a[0] is a HiFi read identifier.
                if int(a[1]) == 0:
                    Read_1 = hifiwrangler0.Read(a[0],i)
                    # add the mapped spatial coordinate/spot identifer,a[2],to this new instance. 
                    Read_1.add_spot(a[2])
                    #
                    # the following codes work to extract
                    # (i)the raw identifier of the mapped spatial coordinate/spot(barcode),a[2]
                    # (ii)the N value - the number of barcodes have the same sequence with a[2]
                    my_re = re.compile(r'(\S+)_(\d+)')
                    my_re_res = my_re.search(a[2])
                    a_raw_spot = my_re_res.group(1)
                    Read_1.add_raw_spot(a_raw_spot)
                    # This HiFi read identifier is stored in a dictionary so this instance can be retrived easily.
                    # Are there any other more efficient way? 
                    HiFi[a[0]] = Read_1
                    N = my_re_res.group(2)
                    # Use a dictionary Spot_dic to store each N for each spot identifier.
                    Spot_dic[a_raw_spot] = N
                else:
                    # There should be a 256 in the second field (flag) in this line.
                    # https://en.wikipedia.org/wiki/SAM_(file_format)#Format 
                    if a[0] in HiFi:
                        as1 = int(i.replace("AS:i:",""))
                        if as1 >= int(HiFi[a[0]].AS):
                            # every HiFi read has flag-256 aligned spot may have flag-0 aligned spot and therefore it has been stored in HiFi as an instance. 
                            # We could add this flag-256 aligned spot to this HiFi read using the add_spot method.
                            HiFi[a[0]].add_spot(a[2])   
                            my_re = re.compile(r'(\S+)_(\d+)')
                            my_re_res = my_re.search(a[2])
                            a_raw_spot = my_re_res.group(1)
                            HiFi[a[0]].add_raw_spot(a_raw_spot)                    
                            N = my_re_res.group(2)
                            Spot_dic[a_raw_spot] = N
                    else:
                        # in such cases, a hifi read was aligned to spots in a forward/forward way WITHOUT a flag of 0 but only flag(s) of 256.
                        # in the SAM file, aligned spots were ranked by the AS field. 
                        Read_1 = hifiwrangler0.Read(a[0],i)
                        Read_1.add_spot(a[2])
                        my_re = re.compile(r'(\S+)_(\d+)')
                        my_re_res = my_re.search(a[2])
                        a_raw_spot = my_re_res.group(1)
                        Read_1.add_raw_spot(a_raw_spot)
                        HiFi[a[0]] = Read_1
                        N = my_re_res.group(2)
                        Spot_dic[a_raw_spot] = N

# the key in Spot_dic has no the '_N' part, they are the raw spot identifier
for Spot in Spot_dic.keys():
    Spot_1 = hifiwrangler0.Spot(Spot)   
    Spot_obj[Spot] = Spot_1

with open(spot_duplic_info) as file:
    for line in file:
        a = line.split("\t")
        if len(a) == 2:
            if a[0] in Spot_dic:
                # This spot a[0] has been aligned by some HiFi READ. 
                # So it could be found in Spot_dic and also, Spot_obj
                # Now we add a[1], which is not used for BWA-MEM mapping and absent in spot_to_hifi_sam but has the SAME sequence as a[0]
                # we need to remove '\n' from a[1]
                a_raw_spot = a[1].strip()
                Spot_obj[a[0]].add_dup_spot(a_raw_spot)

for k in HiFi.keys():
    hifi = HiFi[k]
    N = 0
    for s in set(hifi.spots):
        my_re = re.compile(r'(\S+)_(\d+)')
        my_re_res = my_re.search(s)
        N = N + int(my_re_res.group(2))
    if N <= Nmax or Nmax == 0:
        # raw spots: spot identifier without the '_N' part.
        for s in set(hifi.raw_spots):
            s0 = s
            s0_fields = s0.split(":")
            print(hifi.name,s0_fields[4],s0_fields[6],s0_fields[5],N,hifi.AS,end="\n",sep="\t")
            for s1 in Spot_obj[s0].duplicates:
                s1_fields = s1.split(":")
                print(hifi.name,s1_fields[4],s1_fields[6],s1_fields[5],N,hifi.AS,end="\n",sep="\t")

