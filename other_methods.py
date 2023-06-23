import EoN
import networkx as nx
import TARank
import time
from scipy.stats import kendalltau


def load_txt(path):
    file = open(path)
    Graph = nx.Graph()
    for lines in file:
        line = lines.rstrip("\n")
        node = line.split(" ")
        Graph.add_edge(node[0], node[1])
    mapping = dict(zip(Graph, range(1, nx.number_of_nodes(Graph) + 1)))
    Graph = nx.relabel_nodes(Graph, mapping)
    return Graph


def get_transmission_rate(Graph):
    D = 0
    D_2 = 0
    for each_node in nx.nodes(Graph):
        D = D + nx.degree(Graph, each_node)
        D_2 = D_2 + nx.degree(Graph, each_node) ** 2
    k = D / nx.number_of_nodes(Graph)
    k_2 = D_2 / nx.number_of_nodes(Graph)
    a = 1.5
    b = a * (k / (k_2 - k))
    return b


def Kendall_list(Graph, transmission):
    R_list = list()
    R_dict = dict()
    for each_node in range(1, nx.number_of_nodes(Graph) + 1):
        print("节点", each_node)
        iteration = 1
        R = 0
        while iteration <= 1000:
            t_temp, S_temp, I_temp, R_temp = EoN.basic_discrete_SIR(Graph, transmission, initial_infecteds=each_node)
            R = R + R_temp[len(R_temp) - 1]
            iteration = iteration + 1
        R = R / 1000
        R_list.append(R)
        R_dict[each_node] = R
    result_list = list()
    R_dict = sorted(R_dict.items(), key=lambda x: x[1], reverse=True)
    for each in R_dict:
        result_list.append(each[0])
    return R_list, result_list


# 读取文件
file_path = "C:\\Users\\DELL\\Desktop\\test_centrality\\network\\celegans_Phenotypes"
G = load_txt(file_path + ".txt")
t = time.perf_counter()
result_list, result_rank = TARank.process(G)
print(f'coast:{time.perf_counter() - t:.8f}s')
# 获取SIR模型R列表
# p = get_transmission_rate(G)
# kendall_list, kendall_rank = Kendall_list(G, p)
# tau, p = kendalltau(result_list, kendall_list)
# print(tau)


