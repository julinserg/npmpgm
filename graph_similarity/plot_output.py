import numpy as np
import matplotlib.pyplot as plt
import os

import util_fns as uf
import get_graph_props as gp

def plot_avg_degree_dist(graphs, params, ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    n = params['n']
    r = params['r']
    p = params['p']
    # guess some max_deg to keep
    # length of deg_pdf consistent
    max_deg = n
    avg_degs_pdf = 0
    for graph in graphs:
        degs = gp.get_degrees(graph)
        # max_deg = int(np.max(degs))
        degs_pdf = np.empty(max_deg)
        for i in range(max_deg):
            degs_pdf[i] = degs[degs == i+1].size
        avg_degs_pdf = avg_degs_pdf + degs_pdf
    avg_degs_pdf = avg_degs_pdf/(float(n)*len(graphs))
    ax.scatter(np.arange(1, n+1), avg_degs_pdf)
    ax.set_title("Chung-Lu degree dist., p=" + str(p) + ", r=" + str(r))
    ax.set_ylim(bottom=0)
    ax.set_xlim((0, max_deg))
    plt.savefig("/home/oakridge/holiday/workspace/graph_similarity/figs/chunglu_degdist/dd_p" + str(p) + "_r" + str(r) + ".png")
    print "saved in: ", "/home/oakridge/holiday/workspace/graph_similarity/figs/chunglu_degdist/dd_p" + str(p) + "_r" + str(r) + ".png"
    

def plot_chunglu_deg_dists():
    from subprocess import call
    import os
    r_vals = [0, 0.03, 0.06, 0.09]
    p_vals = [0.7, 0.8, 0.9, 1.0]
    # r_vals = [0, 0.03]
    # p_vals = [0.7]
    for r in r_vals:
        for p in p_vals:
            call([os.getcwd() + "/chunglu_gen", str(p), str(r)])
            graphs = []
            for infile in os.listdir(os.getcwd() + "/output_data"):
                graph, params = uf.get_data(os.getcwd() + "/output_data/" + infile, header_rows=1)
                graphs.append(graph)
            plot_avg_degree_dist(graphs, params)

def plot_chunglu_coloring(eigvals, eigvects, graph_params, output_filename=None, ax=None, nplots=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    eigvals = np.abs(eigvals)
    sorted_indices = np.argsort(eigvals)
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[sorted_indices, :]
    pvals = graph_params[:, 0]
    rvals = graph_params[:, 1]
    # the top eigenvector is basically
    # a vector of ones, do not plot
    # also, must take the "bottom"
    # vectors as the eigenvalues were sorted
    # low to high
    if output_filename is None:
        output_filename = ""
    else:
        output_filename = output_filename + "_"
    ax.set_xlim((np.min(pvals), np.max(pvals)))
    ax.set_ylim((np.min(rvals), np.max(rvals)))
    ax.set_xlabel('p')
    ax.set_ylabel('r')
    if nplots is None:
        nplots = eigvals.shape[0] - 1
    for i in range(2, nplots+2):
        ax.scatter(pvals, rvals, c=eigvects[-i], s=50, lw=0, alpha=0.7)
        ax.set_title("Chung-Lu ensemble coloring: eigenvector " + str(i-1))
        plt.savefig("./figs/chunglu_coloring/" + output_filename + "eigvect_" + str(i-1) + ".png")

def plot_erdosrenyi_coloring(eigvals, eigvects, graph_params, output_filename=None, ax=None, nplots=4):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    eigvals = np.abs(eigvals)
    sorted_indices = np.argsort(eigvals)
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[sorted_indices, :]
    pvals = graph_params
    # the top eigenvector is basically
    # a vector of ones, do not plot
    # also, must take the "bottom"
    # vectors as the eigenvalues were sorted
    # low to high
    if output_filename is None:
        output_filename = ""
    else:
        output_filename = output_filename + "_"
    ax.set_xlim((np.min(pvals), np.max(pvals)))
    ax.set_xlabel('p')
    ax.set_yticks([])
    ax.set_yticklabels("")
    n = pvals.shape[0]
    for i in range(2, nplots+2):
        ax.scatter(pvals, np.zeros(n), c=eigvects[-i], s=50, lw=0, alpha=0.7)
        ax.set_title("Erdos-Renyi ensemble coloring: eigenvector " + str(i-1))
        plt.savefig("./figs/erdosrenyi_coloring/" + output_filename + "eigvect_ " + str(i-1) + ".png")

def plot_erdosrenyi_embedding(eigvals, eigvects, graph_params, output_filename=None, ax=None, maxindex=5, t=0):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    eigvals = np.abs(eigvals)
    sorted_indices = np.argsort(eigvals)
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[sorted_indices, :]
    pvals = graph_params
    # the top eigenvector is basically
    # a vector of ones, do not plot
    # also, must take the "bottom"
    # vectors as the eigenvalues were sorted
    # low to high
    eigvects_to_plot = np.array([eigvects[-i,:] for i in range(2, maxindex+1)])
    eigvals_to_plot = np.array([eigvals[-i] for i in range(2, maxindex+1)])
    nvects = maxindex - 1
    if output_filename is None:
        output_filename = ""
    else:
        output_filename = output_filename + "_"
    for i in range(nvects):
        for j in range(i+1, nvects):
            xvals = np.power(eigvals_to_plot[i], t)*eigvects_to_plot[i,:]
            yvals = np.power(eigvals_to_plot[j], t)*eigvects_to_plot[j,:]
            ax.scatter(xvals , yvals, c=graph_params, lw=0, alpha=0.7)
            ax.set_xlim((np.min(xvals), np.max(xvals)))
            ax.set_ylim((np.min(yvals), np.max(yvals)))
            ax.set_xlabel('Eigvect ' + str(i+1))
            ax.set_ylabel('Eigvect ' + str(j+1))
            ax.set_title('Erdos-Renyi embedding')
            plt.savefig("./figs/erdosrenyi_embedding/" + output_filename + "eigvects_ " + str(i+1) + str(j+1) + ".png")


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='*')
    parser.add_argument('--degree-dist', '-dd', action='store_true', default=False)
    parser.add_argument('--chunglu-deg-dists', '-cdd', action='store_true', default=False)
    parser.add_argument('--chunglu-coloring', '-cc', action='store_true', default=False)
    parser.add_argument('--erdosrenyi-coloring', '-ec', action='store_true', default=False)
    parser.add_argument('--erdosrenyi-embedding', '-ee', action='store_true', default=False)
    parser.add_argument('--output-filename', '-fn', '-ofn', type=str, nargs=1, default="")
    args = parser.parse_args()
    if args.degree_dist:
        graphs = []
        for infile in args.input_files:
            graph, params = uf.get_data(infile, header_rows=1)
            graphs.append(graph)
        plot_avg_degree_dist(graphs, params)
    elif args.chunglu_deg_dists:
        plot_chunglu_deg_dists()
    elif args.chunglu_coloring:
        eigvals = None
        eigvects = None
        graph_params = None
        for infile in args.input_files:
            if "eigvals" in infile:
                eigvals = uf.get_data(infile, header_rows=0)
            if "eigvects" in infile:
                eigvects = uf.get_data(infile, header_rows=0)
            if "graphparams" in infile:
                graph_params = uf.get_data(infile, header_rows=0)
        plot_chunglu_coloring(eigvals, eigvects, graph_params, output_filename=args.output_filename[0])
    elif args.erdosrenyi_embedding:
        eigvals = None
        eigvects = None
        graph_params = None
        for infile in args.input_files:
            if "eigvals" in infile:
                eigvals, params = uf.get_data(infile, header_rows=1)
            if "eigvects" in infile:
                eigvects, params = uf.get_data(infile, header_rows=1)
            if "graphparams" in infile:
                graph_params = uf.get_data(infile, header_rows=0)
        plot_erdosrenyi_embedding(eigvals, eigvects, graph_params, output_filename=args.output_filename[0])
    elif args.erdosrenyi_coloring:
        eigvals = None
        eigvects = None
        graph_params = None
        for infile in args.input_files:
            if "eigvals" in infile:
                eigvals, params = uf.get_data(infile, header_rows=1)
            if "eigvects" in infile:
                eigvects, params = uf.get_data(infile, header_rows=1)
            if "graphparams" in infile:
                graph_params = uf.get_data(infile, header_rows=0)
        plot_erdosrenyi_coloring(eigvals, eigvects, graph_params, output_filename=args.output_filename[0])
