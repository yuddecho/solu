import os

"""
1. 生成 mmseqs 运行脚本，并处理结果: 
    各个数据源去除同一性超过90%的数据
    同一数据源可溶与不可溶数据之间去除同一性为 100% 的数据
    
    去除 tt_solu 中与 pdb_solu 同一性超过 90% 的数据
    
    去除 tt_insolu 中与 pdb_solu 同一性为 100% 的数据
    去除 chang_insolu 中与 pdb_solu 同一性为 100% 的数据
    去除 nesg_insolu 中与 pdb_solu 同一性为 100% 的数据
    
    合并 tt_solu 与 pdb_solu 数据，得到 tt_pdb_solu
    
    去除 tt_pdb_solu 中与 chang_solu 同一性超过 30% 的数据
    去除 tt_insolu 中与 chang_insolu 同一性超过 30% 的数据
    去除 tt_pdb_solu 中与 nesg_solu 同一性超过 30% 的数据
    去除 tt_insolu 中与 nesg_insolu 同一性超过 30% 的数据

"""


def mmseqs_90(file):
    return f'mmseqs easy-cluster ./db/cleaned/{file} ./db/mmseqs/{file[:-6]} ./db/mmseqs/{file[:-6]}_tmp --min-seq-id 0.9'


def mmseqs_search(file_1, file_2, msi, tag=""):
    return f'mmseqs easy-search ./db/mmseqs/{file_1}_rep_seq{tag}.fasta ./db/mmseqs/{file_2}_rep_seq{tag}.fasta ./db/mmseqs/{file_1}_{file_2}.m8 ./db/mmseqs/{file_1}_{file_2}_tmp --min-seq-id {msi}'


def echo(w, msg):
    w.write(f'echo "{msg}" >> log.txt\n')


# data/ == /root/mm/
root_dir = f'../../data'
cleaned_dir = f'{root_dir}/db/cleaned'

mmseqs_script = f'{root_dir}/mmseqs.sh'
mmseqs_script_w = open(mmseqs_script, 'w', encoding='utf-8')

mmseqs_script_w.write('rm -rf ./db/mmseqs/*\n')
mmseqs_script_w.write('mkdir ./db/mmseqs/res\n')
mmseqs_script_w.write('rm -rf ./log.txt\n')
mmseqs_script_w.write('rm ./db/mmseqs/res.tar.gz')
# 1.
# 各个数据源去除同一性超过90%的数据，得到 *_rep_seq.fasta
for item in os.listdir(cleaned_dir):
    mmseqs_script_w.write(f'{mmseqs_90(item)}\n')
    echo(mmseqs_script_w, f'mmseqs 90: {item}')


# 2.
# 同一数据源可溶与不可溶数据之间去除同一性为 100% 的数据
mmseqs_script_w.write(f"{mmseqs_search('tt_solu', 'tt_insolu', 1.0)}\n")
echo(mmseqs_script_w, f'mmseqs 100: tt_solu tt_insolu')

mmseqs_script_w.write(f"{mmseqs_search('chang_solu', 'chang_insolu', 1.0)}\n")
echo(mmseqs_script_w, f'mmseqs 100: chang_solu chang_insolu')

mmseqs_script_w.write(f"{mmseqs_search('nesg_solu', 'nesg_insolu', 1.0)}\n")
echo(mmseqs_script_w, f'mmseqs 100: nesg_solu nesg_insolu')

# 处理结果
mmseqs_script_w.write(
    """cat <<EOF > seq_100.py
import shutil
import os


# 处理 .m8 文件
def get_record(record_file):
    first_name_record = {}
    second_name_record = {}

    with open(record_file, 'r', encoding='utf-8') as r:
        lines = r.readlines()
        for line in lines:
            line = line.strip()
            line_items = line.split('	')
            first_name_record[line_items[0]] = first_name_record.get(line_items[0], 0) + 1
            second_name_record[line_items[1]] = second_name_record.get(line_items[1], 0) + 1
    print(f'{len(first_name_record.keys())}, {len(second_name_record.keys())}')
    return first_name_record, second_name_record


def filter_data(source_file, res_file, record):
    count = 0
    source_total = 0
    with open(res_file, 'w', encoding='utf-8') as w:
        with open(source_file, 'r', encoding='utf-8') as r:
            while True:
                fasta_name = r.readline()
                seq = r.readline()

                fasta_name = fasta_name.strip()
                if not fasta_name:
                    print(f'{source_file}: {source_total}; {res_file}: {count}')
                    break

                source_total += 1

                fasta_name = fasta_name[1:]
                if record.get(fasta_name, 0) != 0:
                    # print(f'{fasta_name}')
                    continue

                w.write(f'>{fasta_name}\\n')
                w.write(seq)

                count += 1


def run(dir, file_name, m8_file, tag):
    source_file = f'{dir}/{file_name}_rep_seq.fasta'
    res = f'{dir}/{file_name}_rep_seq_100.fasta'

    if tag == 0:
        record_m8, _ = get_record(m8_file)
    if tag == 1:
        _, record_m8 = get_record(m8_file)

    filter_data(source_file, res, record_m8)


root = './db/mmseqs'

tt_solu_tt_insolu_m8 = f'{root}/tt_solu_tt_insolu.m8'
chang_solu_chang_insolu_m8 = f'{root}/chang_solu_chang_insolu.m8'
nesg_solu_nesg_insolu_m8 = f'{root}/nesg_solu_nesg_insolu.m8'

args = [
    ['tt_solu', tt_solu_tt_insolu_m8, 0],
    ['tt_insolu', tt_solu_tt_insolu_m8, 1],
    ['chang_solu', chang_solu_chang_insolu_m8, 0],
    ['chang_insolu', chang_solu_chang_insolu_m8, 1],
    ['nesg_solu', nesg_solu_nesg_insolu_m8, 0],
    ['nesg_insolu', nesg_solu_nesg_insolu_m8, 1]
]

for arg in args:
    run(root, arg[0], arg[1], arg[2])

# 复制 ./db/mmseqs/pdb_solu_rep_seq.fasta 到 ./db/mmseqs/pdb_solu_rep_seq_100.fasta
shutil.copy('./db/mmseqs/pdb_solu_rep_seq.fasta', './db/mmseqs/pdb_solu_rep_seq_100.fasta')

EOF

python seq_100.py

"""
)

# *_rep_seq.fasta 经过处理得到 *_rep_seq_100.fasta

# 3.
# 去除 tt_solu 中与 pdb_solu 同一性超过 90% 的数据
# 去除 tt_insolu 中与 pdb_solu 同一性为 100% 的数据
# 去除 chang_insolu 中与 pdb_solu 同一性为 100% 的数据
# 去除 nesg_insolu 中与 pdb_solu 同一性为 100% 的数据
mmseqs_script_w.write(f"{mmseqs_search('tt_solu', 'pdb_solu', 0.9, '_100')}\n")
echo(mmseqs_script_w, f'mmseqs 100 100: tt_solu pdb_solu')

mmseqs_script_w.write(f"{mmseqs_search('tt_insolu', 'pdb_solu', 1.0, '_100')}\n")
echo(mmseqs_script_w, f'mmseqs 100 100: tt_insolu pdb_solu')

mmseqs_script_w.write(f"{mmseqs_search('chang_insolu', 'pdb_solu', 1.0, '_100')}\n")
echo(mmseqs_script_w, f'mmseqs 100 100: chang_insolu pdb_solu')

mmseqs_script_w.write(f"{mmseqs_search('nesg_insolu', 'pdb_solu', 1.0, '_100')}\n")
echo(mmseqs_script_w, f'mmseqs 100 100: nesg_insolu pdb_solu')

# 处理结果
mmseqs_script_w.write(
    """cat <<EOF > seq_100_100.py
import shutil
import os


# 处理 .m8 文件
def get_record(record_file):
    first_name_record = {}
    second_name_record = {}

    with open(record_file, 'r', encoding='utf-8') as r:
        lines = r.readlines()
        for line in lines:
            line = line.strip()
            line_items = line.split('	')
            first_name_record[line_items[0]] = first_name_record.get(line_items[0], 0) + 1
            second_name_record[line_items[1]] = second_name_record.get(line_items[1], 0) + 1
    print(f'{len(first_name_record.keys())}, {len(second_name_record.keys())}')
    return first_name_record, second_name_record


def filter_data(source_file, res_file, record):
    count = 0
    source_total = 0
    with open(res_file, 'w', encoding='utf-8') as w:
        with open(source_file, 'r', encoding='utf-8') as r:
            while True:
                fasta_name = r.readline()
                seq = r.readline()

                fasta_name = fasta_name.strip()
                if not fasta_name:
                    print(f'{source_file}: {source_total}; {res_file}: {count}')
                    break

                source_total += 1

                fasta_name = fasta_name[1:]
                if record.get(fasta_name, 0) != 0:
                    # print(f'{fasta_name}')
                    continue

                w.write(f'>{fasta_name}\\n')
                w.write(seq)

                count += 1
                
         
def run(dir, file_name, m8_file, tag):
    source_file = f'{dir}/{file_name}_rep_seq_100.fasta'
    res = f'{dir}/{file_name}_rep_seq_100_100.fasta'
    
    if tag == 0:
        record_m8, _ = get_record(m8_file)
    if tag == 1:
        _, record_m8 = get_record(m8_file)
    
    filter_data(source_file, res, record_m8)
        

root = './db/mmseqs'

tt_solu_pdb_solu_m8 = f'{root}/tt_solu_pdb_solu.m8'
tt_insolu_pdb_solu_m8 = f'{root}/tt_insolu_pdb_solu.m8'
chang_insolu_pdb_solu_m8 = f'{root}/chang_insolu_pdb_solu.m8'
nesg_insolu_pdb_solu_m8 = f'{root}/nesg_insolu_pdb_solu.m8'

args = [
    ['tt_solu', tt_solu_pdb_solu_m8, 0],
    ['tt_insolu', tt_insolu_pdb_solu_m8, 0],
    ['chang_insolu', chang_insolu_pdb_solu_m8, 0],
    ['nesg_insolu', nesg_insolu_pdb_solu_m8, 0]
]

for arg in args:
    run(root, arg[0], arg[1], arg[2])
    
# 同一格式 _100 -> _100_100
for item_solu in ['pdb_solu', 'chang_solu', 'nesg_solu']:
    shutil.copy(f'./db/mmseqs/{item_solu}_rep_seq_100.fasta', f'./db/mmseqs/{item_solu}_rep_seq_100_100.fasta')

EOF

python seq_100_100.py

wait

""")

# *_rep_seq_100.fasta 经过处理得到 *_rep_seq_100_100.fasta

# 4.
# 去除 tt_solu 中与 chang_solu 同一性超过 30% 的数据
# 去除 tt_insolu 中与 chang_insolu 同一性超过 30% 的数据
# 去除 tt_solu 中与 nesg_solu 同一性超过 30% 的数据
# 去除 tt_insolu 中与 nesg_insolu 同一性超过 30% 的数据
# 去除 pdb_solu 中与 chang_solu 同一性超过 30% 的数据
# 去除 pdb_solu 中与 nesg_solu 同一性超过 30% 的数据
mmseqs_script_w.write(f"{mmseqs_search('tt_solu', 'chang_solu', 0.3, '_100_100')}\n")
echo(mmseqs_script_w, f'mmseqs 100 100 30: tt_solu chang_solu')

mmseqs_script_w.write(f"{mmseqs_search('tt_insolu', 'chang_insolu', 0.3, '_100_100')}\n")
echo(mmseqs_script_w, f'mmseqs 100 100 30: tt_insolu chang_insolu')

mmseqs_script_w.write(f"{mmseqs_search('tt_solu', 'nesg_solu', 0.3, '_100_100')}\n")
echo(mmseqs_script_w, f'mmseqs 100 100 30: tt_solu nesg_solu')

mmseqs_script_w.write(f"{mmseqs_search('tt_insolu', 'nesg_insolu', 0.3, '_100_100')}\n")
echo(mmseqs_script_w, f'mmseqs 100 100 30: tt_insolu nesg_insolu')

mmseqs_script_w.write(f"{mmseqs_search('pdb_solu', 'chang_solu', 0.3, '_100_100')}\n")
echo(mmseqs_script_w, f'mmseqs 100 100 30: pdb_solu chang_solu')

mmseqs_script_w.write(f"{mmseqs_search('pdb_solu', 'nesg_solu', 0.3, '_100_100')}\n")
echo(mmseqs_script_w, f'mmseqs 100 100 30: pdb_solu nesg_solu')

# 处理结果
mmseqs_script_w.write(
    """cat <<EOF > seq_100_100_30.py
import shutil
import os


# 处理 .m8 文件
def get_record(record_file):
    first_name_record = {}
    second_name_record = {}

    with open(record_file, 'r', encoding='utf-8') as r:
        lines = r.readlines()
        for line in lines:
            line = line.strip()
            line_items = line.split('	')
            first_name_record[line_items[0]] = first_name_record.get(line_items[0], 0) + 1
            second_name_record[line_items[1]] = second_name_record.get(line_items[1], 0) + 1
    print(f'{len(first_name_record.keys())}, {len(second_name_record.keys())}')
    return first_name_record, second_name_record


def filter_data(source_file, res_file, record, record2):
    count = 0
    source_total = 0
    with open(res_file, 'w', encoding='utf-8') as w:
        with open(source_file, 'r', encoding='utf-8') as r:
            while True:
                fasta_name = r.readline()
                seq = r.readline()

                fasta_name = fasta_name.strip()
                if not fasta_name:
                    print(f'{source_file}: {source_total}; {res_file}: {count}')
                    break

                source_total += 1

                fasta_name = fasta_name[1:]
                
                # 同时筛选
                if record.get(fasta_name, 0) != 0:
                    # print(f'{fasta_name}')
                    continue
                    
                if record2.get(fasta_name, 0) != 0:
                    # print(f'{fasta_name}')
                    continue

                w.write(f'>{fasta_name}\\n')
                w.write(seq)

                count += 1


def run(dir, file_name, m8_file, m8_file2, tag):
    source_file = f'{dir}/{file_name}_rep_seq_100_100.fasta'
    res = f'{dir}/{file_name}_rep_seq_100_100_30.fasta'

    if tag == 0:
        record_m8, _ = get_record(m8_file)
        record_m82, _ = get_record(m8_file2)
        
    if tag == 1:
        _, record_m8 = get_record(m8_file)
        _, record_m82 = get_record(m8_file2)
        

    filter_data(source_file, res, record_m8, record_m82)


root = './db/mmseqs'

tt_solu_chang_solu_m8 = f'{root}/tt_solu_chang_solu.m8'
tt_insolu_chang_insolu_m8 = f'{root}/tt_insolu_chang_insolu.m8'
tt_solu_nesg_solu_m8 = f'{root}/tt_solu_nesg_solu.m8'
tt_insolu_nesg_insolu_m8 = f'{root}/tt_insolu_nesg_insolu.m8'

pdb_solu_chang_solu_m8 = f'{root}/pdb_solu_chang_solu.m8'
pdb_solu_nesg_solu_m8 = f'{root}/pdb_solu_nesg_solu.m8'

args = [
    ['tt_solu', tt_solu_chang_solu_m8, tt_solu_nesg_solu_m8, 0],
    ['tt_insolu', tt_insolu_chang_insolu_m8, tt_insolu_nesg_insolu_m8, 0],
    ['pdb_solu', pdb_solu_chang_solu_m8, pdb_solu_nesg_solu_m8, 0]
]

for arg in args:
    run(root, arg[0], arg[1], arg[2], arg[3])

EOF

python seq_100_100_30.py

wait

"""
)

# 处理结果
mmseqs_script_w.write(
    """cat <<EOF > res.py
import shutil
import os

root = './db/mmseqs'

# 合并 pdb_solu 和 ttdb_solu
pdb_solu = f'{root}/pdb_solu_rep_seq_100_100_30.fasta'
tt_solu = f'{root}/tt_solu_rep_seq_100_100_30.fasta'
tt_pdb_solu = f'{root}/res/tt_pdb_solu.fasta'

with open(tt_pdb_solu, 'w', encoding='utf-8') as w:
    for file in [pdb_solu, tt_solu]:
        with open(file, 'r', encoding='utf-8') as r:
            while True:
                fasta = r.readline()
                seq = r.readline()
                
                fasta = fasta.strip()
                if not fasta:
                    break
                    
                w.write(f'{fasta}\\n')
                w.write(seq)
            
# 规范名称
test_set = ['chang_solu', 'chang_insolu', 'nesg_solu', 'nesg_insolu']
for item in test_set:
    shutil.copy(f'{root}/{item}_rep_seq_100_100.fasta', f'{root}/res/{item}.fasta')
    
shutil.copy(f'{root}/tt_insolu_rep_seq_100_100_30.fasta', f'{root}/res/tt_insolu.fasta')

EOF

python res.py

wait

"""
)

mmseqs_script_w.write('tar -zcvf res.tar.gz res/ \n')

# 清理文件
mmseqs_script_w.write('rm -rf ./db/mmseqs/*_tmp\n')
mmseqs_script_w.write('rm -rf ./db/mmseqs/*_all_seqs.fasta\n')
mmseqs_script_w.write('rm -rf ./db/mmseqs/*_cluster.tsv\n')

# mmseqs_script_w.write('rm -rf ./db/mmseqs/tt_solu_rep_seq.fasta\n')
# mmseqs_script_w.write('rm -rf ./db/mmseqs/pdb_solu_rep_seq.fasta\n')

mmseqs_script_w.write('rm -rf ./db/mmseqs/*.m8\n')
# mmseqs_script_w.write('rm -rf ./*.py\n')

mmseqs_script_w.close()
