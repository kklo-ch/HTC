import time

import EoN
import networkx as nx
from scipy.stats import kendalltau


def load_txt(path):
    file = open(path)
    Graph = nx.Graph()
    for lines in file:
        line = lines.rstrip("\n")
        node = line.split(" ")
        Graph.add_edge(node[0], node[1])
    mapping = dict(zip(Graph, range(1, nx.number_of_nodes(Graph) + 1)))
    print(mapping)
    Graph = nx.relabel_nodes(Graph, mapping)
    return Graph, mapping


def load_pg_exp(path, node_map):
    file = open(path)
    pg_exp_dict = dict()
    result_list = list()
    result_rank = list()
    for lines in file:
        line = lines.rstrip("\n")
        str_temp = line.split("\t")
        node = node_map[str_temp[0]]
        score = float(str_temp[1])
        pg_exp_dict[node] = score
    pg_exp_dict = sorted(pg_exp_dict.items(), key=lambda x: x[0], reverse=False)
    for each in pg_exp_dict:
        result_list.append(each[1])
    pg_exp_dict = sorted(dict(pg_exp_dict).items(), key=lambda x: x[1], reverse=True)
    for each in pg_exp_dict:
        result_rank.append(each[0])
    return result_list, result_rank


def load_pg_geo(path, node_map):
    file = open(path)
    pg_geo_dict = dict()
    result_list = list()
    result_rank = list()
    for lines in file:
        line = lines.rstrip("\n")
        str_temp = line.split("\t")
        node = node_map[str_temp[0]]
        score = float(str_temp[1])
        pg_geo_dict[node] = score
    pg_geo_dict = sorted(pg_geo_dict.items(), key=lambda x: x[0], reverse=False)
    for each in pg_geo_dict:
        result_list.append(each[1])
    pg_geo_dict = sorted(dict(pg_geo_dict).items(), key=lambda x: x[1], reverse=True)
    for each in pg_geo_dict:
        result_rank.append(each[0])
    return result_list, result_rank


def degree_centrality(Graph):
    dc_list = list()
    dc_rank = list()
    dc_dict = nx.degree_centrality(Graph)
    dc_dict1 = sorted(dc_dict.items(), key=lambda x: x[0], reverse=False)
    for each in dc_dict1:
        dc_list.append(each[1])
    dc_dict2 = sorted(dc_dict.items(), key=lambda x: x[1], reverse=True)
    for each in dc_dict2:
        dc_rank.append(each[0])
    return dc_list, dc_rank


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
file_path = "C:\\Users\\DELL\\Desktop\\test_centrality\\network\\"
exp_path = "C:\\Users\\DELL\\Desktop\\pg_test\\output\\epg\\"
geo_path = "C:\\Users\\DELL\\Desktop\\pg_test\\output\\gpg\\"
file_name = "celegans_Phenotypes"
G, node_index = load_txt(file_path + file_name + ".txt")
# 获取PG排名
t = time.perf_counter()
pg_exp_list, pg_exp_rank = load_pg_exp(exp_path + file_name + ".txt", node_index)
print(f'coast:{time.perf_counter() - t:.8f}s')
t = time.perf_counter()
pg_geo_list, pg_geo_rank = load_pg_geo(geo_path + file_name + ".txt", node_index)
print(f'coast:{time.perf_counter() - t:.8f}s')
# 获取SIR模型R列表
p = get_transmission_rate(G)
kendall_list, kendall_rank = Kendall_list(G, p)
# 计算Kendall相关系数
tau1, p1 = kendalltau(pg_exp_list, kendall_list)
tau2, p2 = kendalltau(pg_geo_list, kendall_list)
print(tau1)
print(tau2)
# 其他算法
dc_list, dc_rank = degree_centrality(G)
dc_tau, dc_p = kendalltau(dc_list, kendall_list)
print(dc_tau)
