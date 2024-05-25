
path=r'G:\SRR21661978\SRR21661978.out-splice-high.sorted.bed'   #路径
offset=1000 #  每多少行生成一个坐标信息，运存不够的电脑这里设置大一点，一般都没问题

#############################  代码区  ################################
out=path[0:-3]+'caqidx'
file = open(path,'r')

#编写索引文件 （place,seek）
#{'chr1': [[0, 0], [1300, 32]], 'chr2': [[0, 52]]}
def file_seek_index_maker(file):
    global offset

    file.seek(0)
    retur={}
    last_seek=0
    skip_time=0
    span_end_dict={}

    while 1 == 1:
        now_line=file.readline()
        now_line_lis=now_line.split('\t')

        if '\n' in now_line_lis:
            now_line_lis.remove('\n')
        print('\r'+str(now_line_lis),end='')
        if now_line_lis == [] or now_line_lis == ['']:
            break

        chr_now=now_line_lis[0]
        place=int(now_line_lis[1])

        if chr_now in retur.keys():
            if int(now_line_lis[2]) > int(span_end_dict[chr_now]):
                span_end_dict[chr_now] = int(now_line_lis[2])

            if skip_time >= offset:
                retur[chr_now] += [[place,last_seek]]
                skip_time = 0
            else:
                skip_time += 1
        else:
            skip_time = 0
            retur[chr_now]=[[0,last_seek]]
            span_end_dict[chr_now] = 0

        last_seek=file.tell()
    retur['chr_span_end']=span_end_dict
    return retur

out_file=open(out,'a')
out_file.write(str(file_seek_index_maker(file)))
out_file.close()


