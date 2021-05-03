import json
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd


def readFileAndProcess(filePath):
    file = open(filePath).readlines()
    # each sequence's identifier(ID number + the species name) as key, sequence as value
    fasta_dict = {}
    for line in file:
        if line.startswith(">"):
            lineTemp = line.replace("|", "")
            name = lineTemp.replace("\n", "")
            # can add in a space or | or whatever in between [0] and [1]
        else:
            fasta_dict[name] = line.replace("\n", "")
    return fasta_dict
# calculate the differences of pariwise sequence


def readFileAndProcessB(filePath):
    file = open(filePath).readlines()
    # each sequence's identifier(ID number + the species name) as key, sequence as value
    nameArr = {}
    for line in file:
        lineTemp = line.replace("|", "")
        lineTemp = lineTemp.replace("\n", "")
        seqArr = lineTemp.split("=")
        for seq in seqArr:
            nameArr[seq] = seqArr
    return nameArr


def compute_difference(seq1, seq2):
    difference = 0
    for a, b in zip(seq1, seq2):
        if a != b:
            difference += 1
    return difference


def getDiversityFactorProcessResults(name_list, seq_list):
    results = []
    length = len(seq_list) - 1
    for i in range(length, 0, -1):
        for j in range(i):
            difference = compute_difference(seq_list[i], seq_list[j])
            diversity_factor = (difference/600) if (difference > 0) else 0
            results.append([name_list[i], name_list[j], diversity_factor])
    return results


def printOkSeq(results, input_diversity_factor):
    meet_conditions_results = []
    meet_conditions_seq_combination = ""
    for result in results:
        if result[2] < input_diversity_factor:
            meet_conditions_seq_combination += result[0] if (
                result[0] not in meet_conditions_seq_combination) else ""
            meet_conditions_seq_combination += result[1] if (
                result[1] not in meet_conditions_seq_combination) else ""
            meet_conditions_results.append(
                [result[0], result[1], str(result[2]*100) + "%"])
    # print("those sequences meet the condition:")
    # print(meet_conditions_results)
    # print("combine those sequences which meet the condition:")
    # print(meet_conditions_seq_combination)
    with open("meet_conditions_results.txt", "w") as out1:
        for m in meet_conditions_results:
            for n in m:
                out1.write(n)
                out1.write(" ")
            out1.write("\n")
        out1.close()

    return meet_conditions_results


def writeFile(filePath, text):
    with open(filePath, 'a') as f:
        f.write('\n' + text)


def processRepeat(meet_conditions_results, nd_separate):
    data = []
    isRepeatArr = {}
    for item in meet_conditions_results:
        # writeFile(filePath, str(item[0])+","+str(item[1])+","+str(item[2]))
        data.append(
            [item[0], item[1], {'weight': float(item[2].replace("%", ""))}])
        if item[0] in nd_separate:
            for seq in nd_separate[item[0]]:
                text = str(seq)+","+str(item[1])+","+str(item[2])
                if text not in isRepeatArr:
                    isRepeatArr[text] = True
                    data.append(
                        [seq, item[1], {'weight': float(item[2].replace("%", ""))}])
                    # writeFile(filePath, text)
                if len(nd_separate[item[0]]) > 1:
                    for seqB in nd_separate[item[0]]:
                        text = str(seq)+","+str(seqB)+","+str('0%')
                        if text not in isRepeatArr:
                            isRepeatArr[text] = True
                            data.append(
                                [seq, seqB, {'weight': 0}])
                            # writeFile(filePath, text)
        elif item[1] in nd_separate:
            for seq in nd_separate[item[1]]:
                text = str(item[0])+","+str(seq)+","+str(item[2])
                if text not in isRepeatArr:
                    isRepeatArr[text] = True
                    data.append(
                        [item[0], seq, {'weight': float(item[2].replace("%", ""))}])
                    # writeFile(filePath, text)
                if len(nd_separate[item[1]]) > 1:
                    for seqB in nd_separate[item[1]]:
                        text = str(seq)+","+str(seqB)+","+str('0%')
                        if text not in isRepeatArr:
                            isRepeatArr[text] = True
                            data.append(
                                [seq, seqB, {'weight': 0}])
                            # writeFile(filePath, text)
    return data


def generatePic(data, picName):
    G = nx.Graph()  # 创建空图，无向图
    G.add_edges_from(data)
    nx.write_gexf(G, picName + '.gexf')
    # plt.figure(3, figsize=(655, 655))  # 这里控制画布的大小，可以说改变整张图的布局
    # plt.subplot(111)
    # pos = nx.spring_layout(G, iterations=20)
    # # nx.draw(G1_LCC, node_color="red", edge_color="grey", node_size="20")
    # nx.draw(G, pos, edge_color="grey", node_size=10)  # 画图，设置节点大小
    # nx.draw_networkx_labels(G, pos,
    #                         font_size=1)  # 将desc属性，显示在节点上
    # edge_labels = nx.get_edge_attributes(G, 'name')  # 获取边的name属性，
    # nx.draw_networkx_edge_labels(
    #     G, pos, edge_labels=edge_labels, font_size=1)  # 将name属性，显示在边上
    # # 检测圆环
    # print(nx.cycle_basis(G.to_undirected()))
    # # plt.savefig(picName+".png")
    # # plt.show()
    # nx.write_gexf(G, "xxxxx")


def main():
    # clustering the fasta sequences
    fasta_dict = readFileAndProcess("end.fasta")
    writeFile("xxx.json", json.dumps(fasta_dict))
    name_list = list(fasta_dict.keys())
    seq_list = list(fasta_dict.values())
    results = getDiversityFactorProcessResults(name_list, seq_list)
    nd_separate = readFileAndProcessB(
        "new.txt")
    while True:
        input_diversity_factor = input("Please enter max diversity factor:")
        input_diversity_factor = float(input_diversity_factor)
        meet_conditions_results = printOkSeq(results, input_diversity_factor)
        data = processRepeat(meet_conditions_results, nd_separate)
        generatePic(data, "pic"+str(input_diversity_factor))


main()
