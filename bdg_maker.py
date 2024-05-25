
file_path_list=[r'G:\SRR21661978\SRR21661978.out-splice-high.sorted.bed',r'G:\SRR21661978\SRR21661978.out-fanse2high.sorted.bed']
out_path=r'G:\try5.25.bdg'         #必填--普通bdg文件输出路径
out_dict=r'G:\try5.25.dict'          #非必填--生成非官方格式的文件，不希望生成请写为out_dict=r''
how_many_legth_once=100000    #每次统计多少长度，统计完一次写入一次
fail_time_max=25     #读取文件读取到几条不符合就跳过了


##########################################  代码区  #########################################################
class bed_file:
    bed_file_path=''
    bed_file=''
    caqidx_path=''
    caqidx=False


for file_path in file_path_list:
    caqidx_path=file_path[0:-3]+'caqidx'
    globals()[str(file_path)+'____file_class']=bed_file()
    globals()[str(file_path) + '____file_class'].bed_file_path=str(file_path)
    globals()[str(file_path) + '____file_class'].caqidx_path=caqidx_path

    try:
        caqfile=open(caqidx_path,'r')
        globals()[str(file_path)+'____file_class'].caqidx=eval(caqfile.read())
    except:
        globals()[str(file_path) + '____file_class'].caqidx=False


end_place_dict={}
for file_path in file_path_list:
    file=open(globals()[str(file_path)+'____file_class'].caqidx_path,'r')
    now_end_dict=eval(file.read())['chr_span_end']  #{'chr1': 24342616, 'chr2': 29647370, 'chr3': 38947415, 'chr4': 34370798, 'chr5': 29512709, 'chr6': 39889241, 'chr7': 23941464}
    for k in now_end_dict.keys():
        if str(k) not in end_place_dict.keys():           #不存在当前染色体数据，进行创建
            end_place_dict[str(k)]=int(now_end_dict[k])
        else:
            if now_end_dict[k] > end_place_dict[str(k)]:     #当前染色体数据进行更新
                end_place_dict[str(k)] = now_end_dict[k]
            else:
                pass

print('各染色体结束位置字典：{}'.format(str(end_place_dict)))

def count_dict_updata(lis):
    global count_dict
    lis_start=int(lis[1])
    lis_end=int(lis[2])
    for i in count_dict.keys():
        if int(i) >= lis_start and int(i) < lis_end:
            count_dict[i] += 1



for chr_name in end_place_dict.keys():          #对每个染色体进行操作
    print('现在对{}染色体进行统计'.format(str(chr_name)))
    chr_end=int(end_place_dict[chr_name])
    range_start_now=0

    while range_start_now < chr_end:
        count_dict={}
        for place_int in range(range_start_now,range_start_now+how_many_legth_once):#创建统计字典
            count_dict[place_int]=0
        print('对{}染色体{}-{}段进行统计中'.format(chr_name,range_start_now,range_start_now+how_many_legth_once))
        for file_path in file_path_list:              #分别对列表中文件进行操作
            file_now=open(file_path,'r')

            if globals()[str(file_path) + '____file_class'].caqidx != False:       #如果有caqidx文件的话
                good_start_index_have = 0
                for idx in globals()[str(file_path) + '____file_class'].caqidx[chr_name]:  # 利用caqidx搜寻合适起点
                    if idx[0] > range_start_now:
                        try:
                            good_start_index = last_index
                            good_start_index_have = 1
                            break
                        except:
                            last_index = idx[1]
                            good_start_index = last_index
                            good_start_index_have = 1
                            break
                    else:
                        last_index = idx[1]

                if good_start_index_have != 1:
                    good_start_index = globals()[str(file_path) + '____file_class'].caqidx[chr_name][-1][1]

            else:          #没有索引文件的话，从头开始
                good_start_index=0


            #起始点已经确定好，准备进行下一步操作

            file_now.seek(good_start_index)

            while 1 == 1:                       #寻找开头
                now_line=file_now.readline()
                now_lis=now_line.split('\t')
                found_need_chr=False
                if len(now_lis) < 5:
                    print('读取停止，当前行列表{}'.format(str(now_lis)))
                    start_success=0
                    break

                if str(now_lis[0]) == str(chr_name) and found_need_chr == False:
                    found_need_chr = True

                if found_need_chr == True and str(now_lis[0]) != str(chr_name):
                    start_success=0
                    break

                if now_lis[0] == str(chr_name) and range_start_now < int(now_lis[2]) and range_start_now+how_many_legth_once > int(now_lis[1]):#发现符合的了
                    count_dict_updata(now_lis)  #对符合的处理函数
                    start_success=1
                    break

            fail_times=0
            stop_sigle=0
            if start_success == 1:
                while stop_sigle != 1:
                    now_line = file_now.readline()
                    now_lis = now_line.split('\t')
                    if now_lis[0] == str(chr_name) and range_start_now < int(now_lis[2]) and range_start_now + how_many_legth_once > int(now_lis[1]): #符合要求
                        count_dict_updata(now_lis)
                        fail_times = 0
                    else:
                        fail_times += 1
                        if fail_times >= fail_time_max:
                            stop_sigle = 1
                            break
                        if now_lis[0] != str(chr_name):
                            stop_sigle = 1
                            break




            file_now.close()
        #统计结束，马上写入
        print('正在写入，现在请不要停止程序')
        out_file=open(out_path,'a')
        for i in count_dict.keys():
            if count_dict[i] != 0 :
                write_line = '{}\t{}\t{}\t{}\n'.format(chr_name,int(i),int(i)+1,str(count_dict[i]))
                out_file.write(write_line)
        out_file.close()


        if out_dict != '':
            out_dict_file = open(out_dict,'a')
            for ii in count_dict.keys():
                if count_dict[ii] != 0:
                    xx = "'{}_{}':{},".format(chr_name, int(ii), count_dict[ii])
                    out_dict_file.write(xx)

            out_dict_file.close()


        print('写入完成')

        range_start_now += how_many_legth_once




