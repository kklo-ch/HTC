import EoN
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import kendalltau


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


def uniform_length(list_array, number):
    if number == 0:
        max_length = 0
        for each_list in list_array:
            if max_length < len(each_list):
                max_length = len(each_list)
        for each_list in list_array:
            times = max_length - len(each_list)
            i = 1
            while i <= times:
                each_list.append(each_list[len(each_list) - 1])
                i = i + 1
        return list_array
    elif number == 1:
        max_length = 0
        new_list = list()
        for each_list in list_array:
            if max_length < len(each_list):
                max_length = len(each_list)
        for each_list in list_array:
            if len(each_list) == max_length:
                new_list = each_list
        return new_list


def SIR_model(Graph, central_nodes, p):
    iteration = 1
    t_sum, S_sum, I_sum, R_sum = EoN.basic_discrete_SIR(Graph, p, initial_infecteds=central_nodes)
    while iteration <= 99:
        t, S, I, R = EoN.basic_discrete_SIR(Graph, p, initial_infecteds=central_nodes)
        if len(I_sum) <= len(I):
            temp1 = list()
            temp2 = list()
            i = 0
            while i < len(I_sum):
                I_temp = I_sum[i] + I[i]
                R_temp = R_sum[i] + R[i]
                temp1.append(I_temp)
                temp2.append(R_temp)
                i = i + 1
            while i < len(I):
                I_temp = I_sum[len(I_sum) - 1] + I[i]
                R_temp = R_sum[len(R_sum) - 1] + R[i]
                temp1.append(I_temp)
                temp2.append(R_temp)
                i = i + 1
            I_sum = temp1
            R_sum = temp2
        elif len(I_sum) > len(I):
            temp1 = list()
            temp2 = list()
            i = 0
            while i < len(I):
                I_temp = I_sum[i] + I[i]
                R_temp = R_sum[i] + R[i]
                temp1.append(I_temp)
                temp2.append(R_temp)
                i = i + 1
            while i < len(I_sum):
                I_temp = I_sum[i] + I[len(I) - 1]
                R_temp = R_sum[i] + R[len(R) - 1]
                temp1.append(I_temp)
                temp2.append(R_temp)
                i = i + 1
            I_sum = temp1
            R_sum = temp2
        iteration = iteration + 1
    Ft_list = list()
    t_list = list()
    j = 0
    while j < len(I_sum):
        Ft_list.append((I_sum[j] + R_sum[j]) / 100)
        t_list.append(j)
        j = j + 1
    return t_list, Ft_list


def SIR_experiment(Graph, p, p_value, dc_list, cc_list, katz, bc_list, ecc_list, pagerank):
    i = 0
    while i < len(p_value):
        plt.clf()
        t1, F1 = SIR_model(Graph, p_value[i], p)
        t2, F2 = SIR_model(Graph, dc_list[i], p)
        t3, F3 = SIR_model(Graph, cc_list[i], p)
        t4, F4 = SIR_model(Graph, katz[i], p)
        t5, F5 = SIR_model(Graph, bc_list[i], p)
        t6, F6 = SIR_model(Graph, ecc_list[i], p)
        t7, F7 = SIR_model(Graph, pagerank[i], p)
        F = list()
        F.append(F1)
        F.append(F2)
        F.append(F3)
        F.append(F4)
        F.append(F5)
        F.append(F6)
        F.append(F7)
        T = list()
        T.append(t1)
        T.append(t2)
        T.append(t3)
        T.append(t4)
        T.append(t5)
        T.append(t6)
        T.append(t7)
        F = uniform_length(F, 0)
        T = uniform_length(T, 1)
        plt.plot(T, F[0], color='red', alpha=0.7, label='p-value', marker='v')
        plt.plot(T, F[1], color='blue', alpha=0.7, label='degree centrality', marker='v')
        plt.plot(T, F[2], color='green', alpha=0.7, label='closeness', marker='v')
        plt.plot(T, F[3], color='orange', alpha=0.7, label='katz', marker='v')
        plt.plot(T, F[4], color='purple', alpha=0.7, label='betweenness', marker='v')
        plt.plot(T, F[5], color='pink', alpha=0.7, label='eigenvector', marker='v')
        plt.plot(T, F[6], color='grey', alpha=0.7, label='pagerank', marker='v')
        plt.xlabel("t")
        plt.ylabel("F(t)")
        plt.title("email_SIR p=" + str(p))
        plt.legend()
        plt.savefig("email_SIR_" + str(i + 1) + ".png", dpi=500, pad_inches=0.0)
        i = i + 1
    return 0


def kendall_experiment(Graph, p, p_value, dc_list, cc_list, katz, bc_list, ecc_list, pagerank):
    R_list = list()
    for each_node in range(1, nx.number_of_nodes(Graph) + 1):
        # print("节点" + str(each_node))
        iteration = 1
        R = 0
        while iteration <= 1000:
            t_temp, S_temp, I_temp, R_temp = EoN.basic_discrete_SIR(Graph, p, initial_infecteds=each_node)
            R = R + R_temp[len(R_temp) - 1] / t_temp[len(t_temp) - 1]
            iteration = iteration + 1
        R = R / 1000
        R_list.append(R)
    # 测试
    print_dict = dict()
    i = 0
    while i < len(R_list):
        print_dict[i + 1] = R_list[i]
        i = i + 1
    print_dict = sorted(print_dict.items(), key=lambda x: x[1], reverse=True)
    print_list = list()
    for each in print_dict:
        print_list.append(each[0])
    print(print_list)
    tau_list = list()
    # 测试
    tau_list.append(abs(kendalltau(p_value, R_list)[0]))
    tau_list.append(abs(kendalltau(dc_list, p_value)[0]))
    tau_list.append(abs(kendalltau(cc_list, p_value)[0]))
    tau_list.append(abs(kendalltau(katz, p_value)[0]))
    tau_list.append(abs(kendalltau(bc_list, p_value)[0]))
    tau_list.append(abs(kendalltau(ecc_list, p_value)[0]))
    tau_list.append(abs(kendalltau(pagerank, p_value)[0]))
    return tau_list
