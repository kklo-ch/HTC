import math
import time

import bigfloat

import EoN
import networkx as nx
from matplotlib import pyplot as plt
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


def draw_graph(Graph):
    pos = nx.spring_layout(Graph)
    nx.draw_networkx_nodes(Graph, pos, node_size=60, label="id")
    nx.draw_networkx_edges(Graph, pos, alpha=0.5)
    nx.draw_networkx_labels(Graph, pos)
    plt.show()


def factorial(n, prec):
    factorial_cache = dict()
    for i in range(n + 1):
        if i == 0:
            factorial_cache[i] = 0
        else:
            factorial_cache[i] = bigfloat.add(factorial_cache[i - 1], bigfloat.log(i, prec), prec)
    return factorial_cache


def quick_comb(n, r, factorial_dict, prec):
    temp = bigfloat.sub(factorial_dict[n], factorial_dict[r], prec)
    return bigfloat.sub(temp, factorial_dict[n - r], prec)


def log_add(a, b, prec):
    # if a > b:
    #     a, b = b, a
    temp1 = bigfloat.sub(a, b, prec)
    temp2 = bigfloat.exp(temp1, prec)
    temp3 = bigfloat.add(1, temp2, prec)
    temp4 = bigfloat.log(temp3, prec)
    return bigfloat.add(b, temp4, prec)


def log_sub_1(b, prec):
    temp1 = bigfloat.exp(b, prec)
    temp2 = bigfloat.sub(1, temp1, prec)
    return bigfloat.log(temp2, prec)


def neighbor_degrees(Graph):
    neighbor_degrees_dict = dict()
    result_list = list()
    for each_node in range(1, nx.number_of_nodes(Graph) + 1):
        neighbor_degree = 0
        for each_neighbor in nx.neighbors(Graph, each_node):
            neighbor_degree = neighbor_degree + nx.degree(Graph, each_neighbor)
        neighbor_degrees_dict[each_node] = neighbor_degree
    sorted_degrees_dict = sorted(neighbor_degrees_dict.items(), key=lambda x: x[1], reverse=True)
    for each in sorted_degrees_dict:
        result_list.append(each[0])
    print(neighbor_degrees_dict)
    return result_list


def p_value_algorithm(Graph):
    # 初始化
    p_value_list = list()
    p_value_dict = dict()
    M = nx.number_of_edges(Graph)
    N = nx.number_of_nodes(Graph)
    neighbor_degrees_dict = dict()
    random_graph_dict = dict()
    factorial_list = factorial(N * (N - 1) // 2)
    T_sum = quick_comb((N - 1) * N // 2, M, factorial_list)
    # 计算每个节点的邻居度之和
    for each_node in range(1, nx.number_of_nodes(Graph) + 1):
        neighbor_degree = 0
        for each_neighbor in nx.neighbors(Graph, each_node):
            neighbor_degree = neighbor_degree + nx.degree(Graph, each_neighbor)
        neighbor_degrees_dict[each_node] = neighbor_degree
    # 将邻居度之和从小到大排序
    sorted_degrees_dict = dict(sorted(neighbor_degrees_dict.items(), key=lambda x: x[1], reverse=False))
    # 初始化邻居度之和为0和1的随机图个数
    random_graph_dict[1] = quick_comb((N - 1) * (N - 2) // 2, M, factorial_list)
    random_graph_dict[2] = math.log(N - 1) + quick_comb((N - 2) * (N - 3) // 2, M - 1, factorial_list)
    random_graph_dict[2] = log_add(random_graph_dict[2], random_graph_dict[1])
    # 完善其他随机图个数
    max_neighbor_degree = max(neighbor_degrees_dict.values())
    for i in range(3, max_neighbor_degree + 1):
        print("度", i)
        edges = i - 1
        if N - 2 >= edges - 1 and M >= edges and (N - 2) * (N - 3) // 2 >= M - edges:
            T1 = math.log(N - 1) + quick_comb(N - 2, edges - 1, factorial_list) + \
                 quick_comb((N - 2) * (N - 3) // 2, M - edges, factorial_list)
        else:
            T1 = 0
        min_num = min(N - 1, edges)
        T2 = 0
        for d in range(2, min_num + 1):
            if M >= edges and d * (d - 1) // 2 + d * (N - 1 - d) >= edges - d and (N - 1 - d) * (N - 2 - d) // 2 >= \
                    M - edges:
                temp = quick_comb(N - 1, d, factorial_list) + \
                       quick_comb(d * (d - 1) // 2 + d * (N - 1 - d), edges - d, factorial_list) + \
                       quick_comb((N - 1 - d) * (N - 2 - d) // 2, M - edges, factorial_list)
            else:
                temp = 0
            if T2 == 0:
                T2 = temp
            elif temp == 0:
                T2 = T2
            else:
                T2 = log_add(T2, temp)
        if T1 == 0:
            T = T2
        elif T2 == 0:
            T = T1
        else:
            T = log_add(T1, T2)
        random_graph_dict[i] = log_add(T, random_graph_dict[i - 1])
    print(random_graph_dict)
    # 从头遍历节点计算p值
    for each_node in range(1, nx.number_of_nodes(Graph) + 1):
        print("节点", each_node)
        p_value = random_graph_dict[neighbor_degrees_dict[each_node]] - bigfloat.BigFloat(T_sum)
        p_value = log_sub_1(p_value)
        p_value_list.append(p_value)
        p_value_dict[each_node] = p_value
        print(bigfloat.exp(p_value))
    # # 建立p值排序
    p_value_dict = sorted(p_value_dict.items(), key=lambda x: x[1], reverse=False)
    result_list = list()
    for each in p_value_dict:
        result_list.append(each[0])
    return p_value_list, result_list


def p_value_optimization(Graph, prec):
    # 初始化
    p_value_list = list()
    p_value_dict = dict()
    M = nx.number_of_edges(Graph)
    N = nx.number_of_nodes(Graph)
    neighbor_degrees_dict = dict()
    random_graph_dict = dict()
    factorial_list = factorial(N * (N - 1) // 2, prec)
    T_sum = quick_comb((N - 1) * N // 2, M, factorial_list, prec)
    # 计算每个节点的邻居度之和
    for each_node in range(1, nx.number_of_nodes(Graph) + 1):
        neighbor_degree = 0
        for each_neighbor in nx.neighbors(Graph, each_node):
            neighbor_degree = neighbor_degree + nx.degree(Graph, each_neighbor)
        neighbor_degrees_dict[each_node] = neighbor_degree
    # 将邻居度之和从小到大排序
    sorted_degrees_dict = dict(sorted(neighbor_degrees_dict.items(), key=lambda x: x[1], reverse=False))
    # 初始化邻居度之和为0和1的随机图个数
    random_graph_dict[1] = quick_comb((N - 1) * (N - 2) // 2, M, factorial_list, prec)
    random_graph_dict_temp1 = bigfloat.log(N - 1, prec)
    random_graph_dict_temp2 = quick_comb((N - 2) * (N - 3) // 2, M - 1, factorial_list, prec)
    random_graph_dict[2] = bigfloat.add(random_graph_dict_temp1, random_graph_dict_temp2, prec)
    random_graph_dict[2] = log_add(random_graph_dict[2], random_graph_dict[1], prec)
    # 完善其他随机图个数
    i = 0
    for each_node in sorted_degrees_dict.keys():
        i = i + 1
        print("节点", i)
        edges = neighbor_degrees_dict[each_node] - 1
        # 两个度之和相差1
        if neighbor_degrees_dict[each_node] - max(random_graph_dict.keys()) == 1:
            if N - 2 >= edges - 1 and M >= edges and (N - 2) * (N - 3) // 2 >= M - edges:
                temp1 = bigfloat.log(N - 1, prec)
                temp2 = quick_comb(N - 2, edges - 1, factorial_list, prec)
                temp3 = quick_comb((N - 2) * (N - 3) // 2, M - edges, factorial_list, prec)
                T1 = bigfloat.add(bigfloat.add(temp1, temp2, prec), temp3, prec)
            else:
                T1 = 0
            min_num = min(N - 1, edges)
            T2 = 0
            for d in range(2, min_num + 1):
                if M >= edges and d * (d - 1) // 2 + d * (N - 1 - d) >= edges - d and \
                        (N - 1 - d) * (N - 2 - d) // 2 >= M - edges:
                    temp1 = quick_comb(N - 1, d, factorial_list, prec)
                    temp2 = quick_comb(d * (d - 1) // 2 + d * (N - 1 - d), edges - d, factorial_list, prec)
                    temp3 = quick_comb((N - 1 - d) * (N - 2 - d) // 2, M - edges, factorial_list, prec)
                    temp = bigfloat.add(bigfloat.add(temp1, temp2, prec), temp3, prec)
                else:
                    temp = 0
                if T2 == 0:
                    T2 = temp
                elif temp == 0:
                    T2 = T2
                else:
                    T2 = log_add(T2, temp, prec)
            if T1 == 0:
                T = T2
            elif T2 == 0:
                T = T1
            else:
                T = log_add(T1, T2, prec)
            random_graph_dict[edges + 1] = log_add(T, random_graph_dict[max(random_graph_dict.keys())], prec)
        # 两个度之和相差2及以上
        elif neighbor_degrees_dict[each_node] - max(random_graph_dict.keys()) >= 2:
            edges1 = max(random_graph_dict.keys())
            edges2 = neighbor_degrees_dict[each_node] - 1
            if N - 2 >= edges1 - 1 and M >= edges1 and (N - 2) * (N - 3) // 2 >= M - edges1:
                temp1 = bigfloat.log(N - 1, prec)
                temp2 = quick_comb(N - 2, edges1 - 1, factorial_list, prec)
                temp3 = quick_comb((N - 2) * (N - 3) // 2, M - edges1, factorial_list, prec)
                T11 = bigfloat.add(bigfloat.add(temp1, temp2, prec), temp3, prec)
            else:
                T11 = 0
            min_num = min(N - 1, edges1)
            T21 = 0
            for d in range(2, min_num + 1):
                if M >= edges1 and d * (d - 1) // 2 + d * (N - 1 - d) >= edges1 - d and (N - 1 - d) * (N - 2 - d) // 2 \
                        >= M - edges1:
                    temp1 = quick_comb(N - 1, d, factorial_list, prec)
                    temp2 = quick_comb(d * (d - 1) // 2 + d * (N - 1 - d), edges1 - d, factorial_list, prec)
                    temp3 = quick_comb((N - 1 - d) * (N - 2 - d) // 2, M - edges1, factorial_list, prec)
                    temp = bigfloat.add(bigfloat.add(temp1, temp2, prec), temp3, prec)
                else:
                    temp = 0
                if T21 == 0:
                    T21 = temp
                elif temp == 0:
                    T21 = T21
                else:
                    T21 = log_add(T21, temp, prec)
            if T11 == 0:
                T1_sum = T21
            elif T21 == 0:
                T1_sum = T11
            else:
                T1_sum = log_add(T11, T21, prec)
            if N - 2 >= edges2 - 1 and M >= edges2 and (N - 2) * (N - 3) // 2 >= M - edges2:
                temp1 = bigfloat.log(N - 1, prec)
                temp2 = quick_comb(N - 2, edges2 - 1, factorial_list, prec)
                temp3 = quick_comb((N - 2) * (N - 3) // 2, M - edges2, factorial_list, prec)
                T12 = bigfloat.add(bigfloat.add(temp1, temp2, prec), temp3, prec)
            else:
                T12 = 0
            min_num = min(N - 1, edges2)
            T22 = 0
            for d in range(2, min_num + 1):
                if M >= edges2 and d * (d - 1) // 2 + d * (N - 1 - d) >= edges2 - d and (N - 1 - d) * (N - 2 - d) // 2 \
                        >= M - edges2:
                    temp1 = quick_comb(N - 1, d, factorial_list, prec)
                    temp2 = quick_comb(d * (d - 1) // 2 + d * (N - 1 - d), edges2 - d, factorial_list, prec)
                    temp3 = quick_comb((N - 1 - d) * (N - 2 - d) // 2, M - edges2, factorial_list, prec)
                    temp = bigfloat.add(bigfloat.add(temp1, temp2, prec), temp3, prec)
                else:
                    temp = 0
                if T22 == 0:
                    T22 = temp
                elif temp == 0:
                    T22 = T22
                else:
                    T22 = log_add(T22, temp, prec)
            if T12 == 0:
                T2_sum = T22
            elif T22 == 0:
                T2_sum = T12
            else:
                T2_sum = log_add(T12, T22, prec)
            number = neighbor_degrees_dict[each_node] - max(random_graph_dict.keys())
            temp_sum = bigfloat.add(min(T1_sum, T2_sum), bigfloat.log(number, prec), prec)
            random_graph_dict[neighbor_degrees_dict[each_node]] = log_add(temp_sum,
                                                                random_graph_dict[max(random_graph_dict.keys())], prec)
    # 从头遍历节点计算p值
    for each_node in range(1, nx.number_of_nodes(Graph) + 1):
        print("节点", each_node)
        p_value = bigfloat.sub(random_graph_dict[neighbor_degrees_dict[each_node]], T_sum, prec)
        # p_value = log_sub_1(p_value, prec)
        p_value_list.append(p_value)
        p_value_dict[each_node] = p_value
        print(bigfloat.exp(p_value, prec))
    # # 建立p值排序
    p_value_dict = sorted(p_value_dict.items(), key=lambda x: x[1], reverse=True)
    result_list = list()
    for each in p_value_dict:
        result_list.append(each[0])
    return p_value_list, result_list


def standard(Graph):
    p_value_list = list()
    p_value_dict = dict()
    N = nx.number_of_nodes(Graph)
    M = nx.number_of_edges(Graph)
    T = math.comb(math.comb(N, 2), M)
    neighbor_degrees_dict = dict()
    random_graph_dict = dict()
    for each_node in range(1, nx.number_of_nodes(Graph) + 1):
        neighbor_degree = 0
        for each_neighbor in nx.neighbors(Graph, each_node):
            neighbor_degree = neighbor_degree + nx.degree(Graph, each_neighbor)
        neighbor_degrees_dict[each_node] = neighbor_degree
    random_graph_dict[1] = math.comb(math.comb(N - 1, 2), M)
    random_graph_dict[2] = (N - 1) * math.comb(math.comb(N - 2, 2), M - 1) + random_graph_dict[1]
    max_neighbor_degree = max(neighbor_degrees_dict.values())
    for i in range(3, max_neighbor_degree + 1):
        edges = i - 1
        if M >= edges:
            T1 = (N - 1) * math.comb(N - 2, edges - 1) * math.comb(math.comb(N - 2, 2), M - edges)
        else:
            T1 = 0
        T2 = 0
        min_num = min(N - 1, edges)
        for d in range(2, min_num + 1):
            if M >= edges:
                temp = math.comb(N - 1, d) * math.comb(math.comb(d, 2) + d * (N - 1 - d), edges - d) * \
                       math.comb(math.comb(N - 1 - d, 2), M - edges)
            else:
                temp = 0
            T2 = T2 + temp
        random_graph_dict[i] = T1 + T2 + random_graph_dict[i - 1]
    for each_node in range(1, nx.number_of_nodes(Graph) + 1):
        print("节点", each_node)
        p_value = decimal.Decimal(1 - decimal.Decimal(random_graph_dict[neighbor_degrees_dict[each_node]] / T))
        p_value_list.append(p_value)
        p_value_dict[each_node] = p_value
        print(p_value)
    p_value_dict = sorted(p_value_dict.items(), key=lambda x: x[1], reverse=False)
    result_list = list()
    for each in p_value_dict:
        result_list.append(each[0])
    return p_value_list, result_list


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


def closeness_centrality(Graph):
    cc_list = list()
    cc_rank = list()
    cc_dict = nx.closeness_centrality(Graph)
    cc_dict1 = sorted(cc_dict.items(), key=lambda x: x[0], reverse=False)
    for each in cc_dict1:
        cc_list.append(each[1])
    cc_dict2 = sorted(cc_dict.items(), key=lambda x: x[1], reverse=True)
    for each in cc_dict2:
        cc_rank.append(each[0])
    return cc_list, cc_rank


def betweenness_centrality(Graph):
    bc_list = list()
    bc_rank = list()
    bc_dict = nx.betweenness_centrality(Graph)
    bc_dict1 = sorted(bc_dict.items(), key=lambda x: x[0], reverse=False)
    for each in bc_dict1:
        bc_list.append(each[1])
    bc_dict2 = sorted(bc_dict.items(), key=lambda x: x[1], reverse=True)
    for each in bc_dict2:
        bc_rank.append(each[0])
    return bc_list, bc_rank


def pagerank(Graph):
    pr_list = list()
    pr_rank = list()
    pr_dict = nx.pagerank(Graph)
    pr_dict1 = sorted(pr_dict.items(), key=lambda x: x[0], reverse=False)
    for each in pr_dict1:
        pr_list.append(each[1])
    pr_dict2 = sorted(pr_dict.items(), key=lambda x: x[1], reverse=True)
    for each in pr_dict2:
        pr_rank.append(each[0])
    return pr_list, pr_rank


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


def expected_num(Graph, distance):
    N = nx.number_of_nodes(Graph)
    M = nx.number_of_edges(Graph)
    p0 = M / (N * (N - 1) // 2)
    temp1 = math.exp(- ((N * p0) ** (distance - 1)) / N)
    temp2 = math.exp(- ((N * p0) ** distance) / N)
    es = (N - 1) * (temp1 - temp2)
    return es


def focus_centrality(Graph):
    fc_list = list()
    fc_dict = dict()
    for each_node in range(1, nx.number_of_nodes(Graph) + 1):
        # 将节点分成n阶邻居
        queue = list()
        queue.append(each_node)
        visited_nodes = list()
        visited_nodes.append(each_node)
        level = 0
        level_num = dict()
        while len(queue) != 0:
            num = len(queue)
            level_num[level] = num
            level = level + 1
            i = 0
            while i < num:
                current_node = queue[0]
                for each_neighbor in nx.neighbors(Graph, current_node):
                    if each_neighbor not in visited_nodes:
                        queue.append(each_neighbor)
                        visited_nodes.append(each_neighbor)
                del(queue[0])
                i = i + 1
        fc = 0
        for each_distance in level_num.keys():
            if each_distance != 0:
                es = expected_num(Graph, each_distance)
                os = level_num[each_distance]
                gs = math.exp(- each_distance ** 2)
                fc = fc + gs * (os - es)
        fc_list.append(fc)
        fc_dict[each_node] = fc
    fc_dict = sorted(fc_dict.items(), key=lambda x: x[1], reverse=True)
    result_list = list()
    for each in fc_dict:
        result_list.append(each[0])
    return fc_list, result_list


if __name__ == '__main__':
    # 读取文件
    # file_path = "C:\\Users\\DELL\\Desktop\\test_centrality\\network\\dolphins"
    # G = load_txt(file_path + ".txt")
    G = nx.read_gml("C:\\Users\\DELL\\Desktop\\test_centrality\\simulation\\n=1000 p=0.01\\1.gml")
    mapping = dict(zip(G, range(1, nx.number_of_nodes(G) + 1)))
    G = nx.relabel_nodes(G, mapping)
    # 画图
    # draw_graph(G)
    # 设置精度
    precision = bigfloat.precision(1000)
    # 获取邻居度之和排名
    # neighbor_rank = neighbor_degrees(G)
    # 计算p值
    print("开始")
    t = time.perf_counter()
    test_list, test_rank = p_value_optimization(G, precision)
    print(f'coast:{time.perf_counter() - t:.8f}s')
    # test_list_1, test_rank_1 = standard(G)
    # 获取SIR模型R列表
    # p = get_transmission_rate(G)
    # kendall_list, kendall_rank = Kendall_list(G, p)
    # 计算Kendall相关系数
    # tau, p = kendalltau(test_list, kendall_list)
    # 显示排名
    print(test_rank)
    # print(test_rank_1)
    # print(neighbor_rank)
    # print(kendall_list)
    # print(tau)
    # 其他算法
    # t = time.perf_counter()
    # dc_list, dc_rank = degree_centrality(G)
    # print(f'coast:{time.perf_counter() - t:.8f}s')
    # t = time.perf_counter()
    # cc_list, cc_rank = closeness_centrality(G)
    # print(f'coast:{time.perf_counter() - t:.8f}s')
    # t = time.perf_counter()
    # bc_list, bc_rank = betweenness_centrality(G)
    # print(f'coast:{time.perf_counter() - t:.8f}s')
    # t = time.perf_counter()
    # pagerank_list, pagerank_rank = pagerank(G)
    # print(f'coast:{time.perf_counter() - t:.8f}s')
    # dc_tau, dc_p = kendalltau(dc_list, kendall_list)
    # cc_tau, cc_p = kendalltau(cc_list, kendall_list)
    # bc_tau, bc_p = kendalltau(bc_list, kendall_list)
    # pagerank_tau, pagerank_p = kendalltau(pagerank_list, kendall_list)
    # print(dc_tau)
    # print(cc_tau)
    # print(bc_tau)
    # print(pagerank_tau)
    # # 对比算法
    # t = time.perf_counter()
    # fc_list, fc_rank = focus_centrality(G)
    # print(f'coast:{time.perf_counter() - t:.8f}s')
    # fc_tau, fc_p = kendalltau(fc_list, kendall_list)
    # print(fc_tau)
