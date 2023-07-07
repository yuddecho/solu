import os
from tqdm import tqdm
import xml.etree.ElementTree as ET


# 1.chang dataset
def chang_cleaning(source, res, tag):
    count = 0
    with open(res, 'w', encoding='utf-8') as w:
        with open(source, 'r', encoding='utf-8') as r:
            last_fasta_name = ''
            while True:
                fasta_name = r.readline()

                fasta_name = fasta_name.strip()
                if not fasta_name:
                    print(f'{res}: {count}')
                    break

                if fasta_name[0] != '>':
                    seq = f'{fasta_name}\n'
                    fasta_name = last_fasta_name
                else:
                    last_fasta_name = fasta_name
                    seq = r.readline()

                # 除去His标签, 除去 两个及以上未知氨基酸, 和 非法字符
                seq_str = ['HHHHHH', 'XX', '*']
                look_dog = False
                for item in seq_str:
                    if item in seq:
                        look_dog = True
                        break

                if look_dog:
                    # print(fasta_name)
                    continue

                count += 1

                if len(fasta_name) > 41:
                    fasta_name = fasta_name[:41]

                w.write(f'{fasta_name}_{tag}_chang{count}\n')
                w.write(f'{seq}')


def run_chang(root):
    # chang
    chang = f'{root}/db'
    chang_solu_csv = f'{chang}/source/chang_solu.csv'
    chang_insolu_csv = f'{chang}/source/chang_insolu.csv'

    chang_solu, chang_insolu = f'{chang}/cleaned/chang_solu.fasta', f'{chang}/cleaned/chang_insolu.fasta'

    chang_cleaning(chang_solu_csv, chang_solu, 'solu')
    chang_cleaning(chang_insolu_csv, chang_insolu, 'insolu')


# 2.nesg
def nesg_cleaning(nesg_csv, nesg_fasta, nesg_solu, nesg_insolu):
    # fasta_id: solu
    nesg_score = {}

    with open(nesg_csv, 'r', encoding='utf-8') as r:
        # titile
        _ = r.readline()

        # data: id,exp,sol
        lines = r.readlines()
        for line in lines:
            line = line.strip()
            line_item = line.split(',')

            fasta_id = line_item[0]
            sol = int(line_item[2])

            res = 1
            if sol == 0:
                res = 0

            nesg_score[fasta_id] = res

    nesg_solu_w = open(nesg_solu, 'w', encoding='utf-8')
    nesg_insolu_w = open(nesg_insolu, 'w', encoding='utf-8')

    with open(nesg_fasta, 'r', encoding='utf-8') as r:
        solu_count, insolu_count = 0, 0
        total_count = 0
        while True:
            fasta = r.readline()
            seq = r.readline()

            fasta = fasta.strip()
            if not fasta:
                print(f'{nesg_solu}: {solu_count}')
                print(f'{nesg_insolu}: {insolu_count}')
                print(f'{nesg_fasta}: {total_count}')
                break

            if fasta[0] != '>':
                print(nesg_solu, fasta)
                raise

            total_count += 1

            # 筛选序列
            seq_str = ['HHHHHH', 'XX']
            look_dog = False
            for item in seq_str:
                if item in seq:
                    look_dog = True
                    break

            if look_dog:
                continue

            # 记录
            fasta = fasta[1:]

            sol = nesg_score.get(fasta, -1)
            if sol == -1:
                print(f'sol = -1: {fasta}')
                continue

            if sol == 1:
                solu_count += 1
                nesg_solu_w.write(f'>{fasta}_solu_nesg{solu_count}\n')
                nesg_solu_w.write(seq)

            if sol == 0:
                insolu_count += 1
                nesg_insolu_w.write(f'>{fasta}_insolu_nesg{insolu_count}\n')
                nesg_insolu_w.write(seq)

    nesg_solu_w.close()
    nesg_insolu_w.close()


def run_nesg(root):
    nesg_csv = f'{root}/db/source/nesg.csv'
    nesg_fasta = f'{root}/db/source/nesg.fasta'

    nesg_solu = f'{root}/db/cleaned/nesg_solu.fasta'
    nesg_insolu = f'{root}/db/cleaned/nesg_insolu.fasta'

    nesg_cleaning(nesg_csv, nesg_fasta, nesg_solu, nesg_insolu)


# 3.pdb
def pdb_cleaning(pdb, pdb_solu):
    # to do
    pdb_solu_w = open(pdb_solu, 'a', encoding='utf-8')

    pdb_solu_count = 0

    # 处理文件
    with open(pdb, 'r', encoding='utf-8') as r:
        while True:
            fasta = r.readline()
            seq = r.readline()

            if not seq:
                break

            try:
                # >2i9z_A mol:protein length:105  Putative septation protein spoVG
                # SNAMKVTDVRLRKIQTDGRMKALVSITLDEAFVIHDLRVIEGNSGLFVAMPSKRTPDGEFRDIAHPINSDMRQEIQDAVMKVYDETDEVIPDKNATSDNEESDEA

                # 去除首尾空格，>
                fasta = fasta.strip()

                if fasta[0] != '>':
                    print(pdb, fasta)
                    raise

                fasta = fasta[1:]

                space_index_1 = fasta.find(' ')
                space_index_2 = fasta.find(' ', space_index_1 + 1)

                # 筛选 蛋白序列
                fasta_name = fasta[:space_index_1]
                seq_type = fasta[space_index_1 + 1: space_index_2]
                if 'protein' not in seq_type:
                    continue

                # 过滤 膜蛋白
                detail = fasta[space_index_2 + 1:]
                if 'membrane protein' in detail:
                    continue

                # ## 规范序列
                seq = seq.strip()

                acids = set(seq)
                if len(acids) <= 4:
                    continue

                # 除去His标签, 除去 两个及以上未知氨基酸, 和 非法字符
                seq_str = ['HHHHHH', 'XX', '*']
                look_dog = False
                for item in seq_str:
                    if item in seq:
                        look_dog = True
                        break

                if look_dog:
                    continue

                pdb_solu_count += 1
                pdb_solu_w.write(f'>{fasta_name}_solu_pdb{pdb_solu_count}\n')
                pdb_solu_w.write(f'{seq}\n')

            except Exception as e:
                print(f'{fasta_name}: e')

    pdb_solu_w.close()

    print(f'{pdb_solu}: {pdb_solu_count}')


def run_pdb(root):
    pdb = f'{root}/db/source/pdb_seqres.txt'
    pdb_solu = f'{root}/db/cleaned/pdb_solu.fasta'

    pdb_cleaning(pdb, pdb_solu)


# ttdb
def ttdb_collecting(tt_xml, tt):
    """
    :param tt_xml: 原始xml文件
    :param tt: data/db/tt.csv'
    :return:
    :note: 由于文件大小5G+，采用lxml的etree.iterpaees()进行解析, 获取原始数据：targetId, status, sequence, annotation
        - targetId: target/targetId
        - status: target/status
        - statusHistory: target/trialList/trial/statusHistoryList/statusHistory/status
        - sequence: target/targetSequenceList/targetSequence/oneLetterCode
        - sequenceType: target/targetSequenceList/targetSequence/sequenceType
        - annotation: 判别是否为 跨膜蛋白, key: transmembrane protein, Claudin(Claudins are a family of nearly two dozen transmembrane proteins)
            - target/targetName
            - target/targetAnnotation
            - target/targetRationale
            - target/remark
            - target/laboratoryList/lab
            - target/targetPartnershipList/targetPartnership/partnershipName
            - target/targetSequenceList/targetSequence/sequenceName
    """
    # 读取 xml 文件 to do
    context = ET.iterparse(tt_xml, events=("start", "end"))

    if os.path.exists(tt):
        os.remove(tt)

    tt_w = open(tt, 'a', encoding='utf-8')

    # 处理 标签合集
    tag_range = ['targetId', 'status', 'targetSequenceList', 'targetName', 'targetAnnotation', 'targetRationale',
                 'remark', 'laboratoryList', 'targetPartnershipList', 'trialList']

    # 处理 从标签里得到的字符串
    def format_label(strs: str) -> str:
        strs = strs.strip()
        strs = strs.replace('\n', '')
        strs = strs.replace('\r', '')
        strs = strs.replace('\t', '')
        strs = strs.replace('  ', ' ')
        strs = strs.replace('  ', ' ')
        strs = strs.replace(',', '.')
        return strs

    # 迭代处理标签
    bar = tqdm(total=130760 + 360433)
    root_tag = None
    for index, (event, elem) in enumerate(context):
        # 保留根标签，用于清除空子标签
        if index == 0:
            root_tag = elem

        # 处理 end 事件遇到的不同标签
        if event == "end":
            if elem.tag == 'target':
                # out: targetId,annotation,targetStatus,targetSeqType:targetseq;...,trialStatus:trialSeqType:trialSeq;...
                targetId, annotation, targetStatus = '', '', ''

                # 分别存放 直接记录组 和 实验记录组
                targetList, trialList = [], []

                for child in elem:
                    if child.tag not in tag_range:
                        continue

                    # target id
                    if child.tag == 'targetId':
                        try:
                            targetId = child.text.strip()
                        except Exception as e:
                            print(f'{targetId}: {e}')
                        finally:
                            continue

                    # target status
                    if child.tag == 'status':
                        try:
                            targetStatus = child.text.strip()
                        except Exception as e:
                            print(f'{targetId}, targetStatus: {e}')
                        finally:
                            continue

                    # sequence target/targetSequenceList/targetSequence/oneLetterCode and sequenceName
                    if child.tag == 'targetSequenceList':
                        for targetSequence in child:
                            if targetSequence.tag == 'targetSequence':
                                # targetSeq:targetStatus
                                targetList_item = []

                                for oneLetterCode in targetSequence:
                                    if oneLetterCode.tag == "oneLetterCode":
                                        try:
                                            text = oneLetterCode.text
                                            seq = "" if text is None else format_label(text)
                                            targetList_item.append(seq)
                                        except Exception as e:
                                            print(f'{targetId}, sequence: {e}')
                                        finally:
                                            continue

                                    # 依然使用 oneLetterCode 作为迭代变量
                                    if oneLetterCode.tag == 'sequenceChemicalType':
                                        try:
                                            text = oneLetterCode.text
                                            seqType = "" if text is None else format_label(text)
                                            targetList_item.append(seqType)
                                        except Exception as e:
                                            print(f'{targetId}, sequenceType: {e}')
                                        finally:
                                            continue

                                    if oneLetterCode.tag == 'sequenceName':
                                        try:
                                            text = oneLetterCode.text
                                            if text is not None:
                                                annotation += f'sequenceName: {format_label(text)}; '
                                        except Exception as e:
                                            print(f'{targetId}, sequenceName: {e}')
                                        finally:
                                            continue

                                targetList.append(targetList_item)
                        continue

                    # trial data
                    # statusHistory: target/trialList/trial/statusHistoryList/statusHistory/status
                    if child.tag == 'trialList':
                        for _trial in child:
                            if _trial.tag == 'trial':
                                # trialStatus:trialSeq:trialSeqType
                                trialList_item = []

                                for _statusHistoryList in _trial:
                                    if _statusHistoryList.tag == 'statusHistoryList':
                                        trialStatus = ''
                                        for _statusHistory in _statusHistoryList:
                                            if _statusHistory.tag == "statusHistory":
                                                for _status in _statusHistory:
                                                    if _status.tag == 'status':
                                                        try:
                                                            text = _status.text
                                                            trialStatus += "" if text is None else f'{format_label(text)};'
                                                        except Exception as e:
                                                            print(f'{targetId}, statusHistory: {e}')
                                                        finally:
                                                            break
                                        trialList_item.append(trialStatus)
                                        continue

                                    if _statusHistoryList.tag == 'trialSequenceList':
                                        # seq:type;...
                                        trialSeq_info = ''

                                        for _trialSequence in _statusHistoryList:
                                            if _trialSequence.tag == 'trialSequence':
                                                for trialSequence_child in _trialSequence:
                                                    if trialSequence_child.tag == 'oneLetterCode':
                                                        try:
                                                            text = trialSequence_child.text
                                                            if text is not None:
                                                                _seq = f'{format_label(text)}'
                                                                trialSeq_info += f'{_seq}+;'
                                                        except Exception as e:
                                                            print(f'{targetId}, trialSequence: {e}')
                                                        finally:
                                                            continue

                                                    if trialSequence_child.tag == 'sequenceChemicalType':
                                                        try:
                                                            text = trialSequence_child.text
                                                            if text is not None:
                                                                _type = f'{format_label(text)}'
                                                                trialSeq_info = trialSeq_info[:-1] + f'{_type};'
                                                        except Exception as e:
                                                            print(f'{targetId}, trialSequenceType: {e}')
                                                        finally:
                                                            continue

                                        trialList_item.append(trialSeq_info)
                                        continue

                                trialList.append(trialList_item)
                                continue
                        continue

                    # annotation
                    if child.tag in ['targetName', 'targetAnnotation', 'targetRationale', 'remark']:
                        try:
                            text = child.text
                            if text is not None:
                                annotation += f'{child.tag}: {format_label(text)}; '
                        except Exception as e:
                            print(f'{targetId}, {child.tag}: {e}')
                        finally:
                            continue

                    # target/laboratoryList/lab
                    if child.tag == 'laboratoryList':
                        for lab in child:
                            if lab.tag == 'lab':
                                try:
                                    text = lab.text
                                    if text is not None:
                                        annotation += f'laboratory: {format_label(text)}; '
                                except Exception as e:
                                    print(f'{targetId}, laboratory: {e}')
                                finally:
                                    continue

                    # target/targetPartnershipList/targetPartnership/partnershipName
                    if child.tag == 'targetPartnershipList':
                        for targetPartnershipList in child:
                            if targetPartnershipList.tag == 'targetPartnership':
                                for targetPartnership in targetPartnershipList:
                                    if targetPartnership.tag == "partnershipName":
                                        try:
                                            text = targetPartnership.text
                                            if text is not None:
                                                annotation += f'{targetPartnership.tag}; '
                                            seq = "" if text is None else format_label(text)
                                        except Exception as e:
                                            print(f'{targetId}, {targetPartnership.tag}: {e}')
                                        finally:
                                            break
                                break
                        continue

                # 清除当前标签内容
                elem.clear()

                # 写入 tt.csv: targetId, annotation, statuse, seq, seqType
                for i in range(len(targetList)):
                    item = targetList[i]
                    try:
                        tt_w.write(f'{targetId},{annotation},{targetStatus},{item[0]},{item[1]}\n')
                        bar.update(1)
                    except Exception as e:
                        print(targetId, item, e)

                for i in range(len(trialList)):
                    item = trialList[i]
                    try:
                        if len(item) == 2:
                            seq = item[1]
                            seq_list = seq.split(';')
                            for seq_item in seq_list:
                                if seq_item != '':
                                    seq_item_list = seq_item.split('+')
                                    tt_w.write(
                                        f'{targetId},{annotation},{item[0]},{seq_item_list[0]},{seq_item_list[1]}\n')
                                    bar.update(1)
                        else:
                            print(targetId, item)
                            raise

                    except Exception as e:
                        print(targetId, item, e)

        # 清除空的子标签
        root_tag.clear()

    tt_w.close()
    bar.close()


def ttdb_cleaning(tt, tt_solu, tt_insolu):
    """
    :param tt:
    :param tt_solu:
    :param tt_insolu:
    :return:
    """
    if os.path.exists(tt_solu):
        os.remove(tt_solu)
        os.remove(tt_insolu)

    # 结果数据
    tt_solu_w = open(tt_solu, 'a', encoding='utf-8')
    tt_insolu_w = open(tt_insolu, 'a', encoding='utf-8')

    tt_solu_count, tt_insolu_count = 0, 0

    # 推断溶解度 key
    status_solu = ['in PDB', 'soluble', 'diffraction', 'crystal structure', 'crystallized',
                   'diffraction-quality crystals', 'NMR assigned', 'native diffraction-data',
                   'phasing diffraction-data', 'HSQC satisfactory', 'NMR structure', 'in BMRB']

    # 逐条处理数据
    with open(tt, 'r', encoding='utf-8') as r:
        lines = r.readlines()
        for line in tqdm(lines):
            line = line.strip()

            # targetId, annotation, status, sequence, sequenceType
            line_list = line.split(',')
            targetId = line_list[0]

            # ## annotation: 去除膜蛋白
            annotation = line_list[1]
            if annotation != '':
                annotation = annotation.lower()
                if 'membrane protein' in annotation or 'claudin' in annotation:
                    continue

            # ## 推断溶解度, 并写入文件
            status = line_list[2]
            status = status.strip()
            solubility = 'insolu'

            if ';' not in status:
                if status == 'work stopped':
                    continue

                if status in status_solu:
                    solubility = 'solu'

            else:
                status = status[:-1]
                status_list = status.split(';')

                for item in status_list:
                    if item in status_solu:
                        solubility = 'solu'
                        break

            # ## 规范 sequence
            sequence = line_list[3]
            sequence = sequence.strip()

            sequenceType = line_list[4]

            # 确保是 氨基酸序列
            if sequenceType != 'protein':
                continue

            acids = set(sequence)
            if len(acids) <= 4:
                continue

            if sequence[0] == '>':
                print(f'{targetId}: {sequence}')
                continue

            # 除去His标签, 除去 两个及以上未知氨基酸, 和 非法字符
            seq_str = ['HHHHHH', 'XX', '_', '.', '*']
            look_dog = False
            for item in seq_str:
                if item in sequence:
                    look_dog = True
                    break

            if look_dog:
                continue

            if solubility == 'solu':
                tt_solu_count += 1
                tt_solu_w.write(f'>{targetId}_{solubility}_tt{tt_solu_count}\n')
                tt_solu_w.write(f'{sequence}\n')
            else:
                tt_insolu_count += 1
                tt_insolu_w.write(f'>{targetId}_{solubility}_tt{tt_insolu_count}\n')
                tt_insolu_w.write(f'{sequence}\n')

    tt_solu_w.close()
    tt_insolu_w.close()

    print(f'TTDB solu: {tt_solu_count}, insolu: {tt_insolu_count}')


def run_ttdb(root):
    # 1. 数据采集
    # tt.xml 包含 TargetTrack DB 所有数据
    # './TargetTrack/TargetTrack-1Jul2017/TargetTrack XML files/tt.xml'
    tt_xml = f'{root}/db/source/tt.xml'

    # 从 tt.xml 文件中获取对应蛋白数据, 解析后存为 tt.csv
    tt = f'{root}/db/source/tt.csv'

    # 采集时间久，执行过一次就不再执行
    if not os.path.exists(tt):
        ttdb_collecting(tt_xml, tt)
    else:
        print(f'collected: {tt}')

    # 2. 数据清洗
    tt_solu, tt_insolu = f'{root}/db/cleaned/tt_solu.fasta', f'{root}/db/cleaned/tt_insolu.fasta'

    if not (os.path.exists(tt_solu) and os.path.exists(tt_insolu)):
        ttdb_cleaning(tt, tt_solu, tt_insolu)
    else:
        print(f'cleaned: {tt_solu}, {tt_insolu}')


def run():
    root = f'../../data'

    # run_chang(root)
    # run_nesg(root)
    # run_pdb(root)
    run_ttdb(root)


if __name__ == '__main__':
    run()
