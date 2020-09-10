# created by Arved Enders-Seidlitz on 02.09.2020
#
# Post processing of Elmer simulation:
# Evaluation of heat fluxes.

import os
import yaml
import meshio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from elements import Node, Line1st, Line2nd, Triangle1st, Triangle2nd


def heat_flux(sim_path, plot=False, save=True, normal_proj=False):
    # import data
    files = os.listdir(sim_path)
    vtu = []
    for f in files:
        ending = f.split('.')[-1]
        if ending == 'msh':
            msh = f
        if ending == 'vtu':
            vtu.append(f)
    msh = meshio.read(sim_path + '/' + msh)
    vtu = meshio.read(sim_path + '/' + vtu[-1])
    with open (sim_path + '/post_processing.yml', 'r') as f:
        data = yaml.safe_load(f)

    # detect interpolation order
    if 'triangle' in msh.cells_dict:
        order = 1
    elif 'triangle6' in msh.cells_dict:
        order = 2

    # extract data from msh and vtu
    if order == 1:
        boundary_ids = msh.cell_data_dict['gmsh:physical']['line']
        nodes = msh.cells_dict['line']
        df_lines = pd.DataFrame({'BoundaryID': boundary_ids, 'Node1': nodes[:, 0], 'Node2': nodes[:, 1]})

        surface_ids = msh.cell_data_dict['gmsh:physical']['triangle']
        nodes = msh.cells_dict['triangle']
        df_elements = pd.DataFrame({'SurfaceID': surface_ids, 'Node1': nodes[:, 0], 'Node2': nodes[:, 1], 'Node3': nodes[:, 2]})

    elif order == 2:
        boundary_ids = msh.cell_data_dict['gmsh:physical']['line3']
        nodes = msh.cells_dict['line3']
        df_lines = pd.DataFrame({'BoundaryID': boundary_ids, 'Node1': nodes[:, 0], 'Node2': nodes[:, 1], 'Node3': nodes[:, 2]})

        surface_ids = msh.cell_data_dict['gmsh:physical']['triangle6']
        nodes = msh.cells_dict['triangle6']
        df_elements = pd.DataFrame({'SurfaceID': surface_ids,
            'Node1': nodes[:, 0], 'Node2': nodes[:, 1], 'Node3': nodes[:, 2],
            'Node4': nodes[:, 3], 'Node5': nodes[:, 4], 'Node6': nodes[:, 5]
            })

    points = vtu.points  # take coordinates from vtu because mesh may have been deformed by PhaseChangeSolver
    temperature = vtu.point_data['temperature']
    df_nodes = pd.DataFrame({'x': points[:, 0], 'y': points[:, 1], 'z': points[:, 2], 'T': temperature[:, 0]})

    # evaluate boundaries
    x_quiver = []  # coordinates for quiver plot
    q_quiver = []  # heat fluxes for quiver plot
    fluxes = {}
    for boundary in data:
        df_l = df_lines.loc[df_lines['BoundaryID'] == data[boundary]['ID']]
        df_e = df_elements.loc[df_elements['SurfaceID'].isin(data[boundary]['BodyIDs'])]
        lmbd = data[boundary]['lambda']
        # find elements on boundary & in body
        element_ids = []
        for _, line in df_l.iterrows():
            if order == 1:
                nodes = [line['Node1'], line['Node2']]
                idx, = df_e.loc[(df_e['Node1'].isin(nodes)) & (df_e['Node2'].isin(nodes)) |
                                (df_e['Node1'].isin(nodes)) & (df_e['Node3'].isin(nodes)) |
                                (df_e['Node2'].isin(nodes)) & (df_e['Node3'].isin(nodes))
                                ].index
            elif order == 2:
                nodes = [line['Node1'], line['Node2'], line['Node3']]
                idx, = df_e.loc[(df_e['Node1'].isin(nodes)) & (df_e['Node4'].isin(nodes)) & (df_e['Node2'].isin(nodes)) |
                                (df_e['Node2'].isin(nodes)) & (df_e['Node5'].isin(nodes)) & (df_e['Node3'].isin(nodes)) |
                                (df_e['Node3'].isin(nodes)) & (df_e['Node6'].isin(nodes)) & (df_e['Node1'].isin(nodes))
                                ].index
            element_ids.append(idx)      
        df_l.insert(3, 'ElementID', element_ids)

        Q = 0
        for _, line in df_l.iterrows():
            element = df_e.loc[line['ElementID'] == df_e.index].iloc[0]
            node = df_nodes.iloc[element['Node1']]
            n1 = Node(node['x'], node['y'], node['z'], node['T'])
            node = df_nodes.iloc[element['Node2']]
            n2 = Node(node['x'], node['y'], node['z'], node['T'])
            node = df_nodes.iloc[element['Node3']]
            n3 = Node(node['x'], node['y'], node['z'], node['T'])
            line_node = df_nodes.iloc[line['Node1']]
            line_n1 = Node(line_node['x'], line_node['y'], line_node['z'])
            line_node = df_nodes.iloc[line['Node2']]
            line_n2 = Node(line_node['x'], line_node['y'], line_node['z'])
            if order == 1:
                t = Triangle1st([n1, n2, n3])
                l = Line1st([line_n1, line_n2])
            if order == 2:
                node = df_nodes.iloc[element['Node4']]
                n4 = Node(node['x'], node['y'], node['z'], node['T'])
                node = df_nodes.iloc[element['Node5']]
                n5 = Node(node['x'], node['y'], node['z'], node['T'])
                node = df_nodes.iloc[element['Node6']]
                n6 = Node(node['x'], node['y'], node['z'], node['T'])
                line_node = df_nodes.iloc[line['Node3']]
                line_n3 = Node(line_node['x'], line_node['y'], line_node['z'])
                t = Triangle2nd([n1, n2, n3, n4, n5, n6])
                l = Line2nd([line_n1, line_n2, line_n3])

            # check orientation of normal, correct if necessary
            if not ((n1.X == l.n1.X).all() or (n1.X == l.n2.X).all()):
                n = n1
            elif not ((n2.X == l.n1.X).all() or (n2.X == l.n2.X).all()):
                n = n2
            elif not ((n3.X == l.n1.X).all() or (n3.X == l.n2.X).all()):
                n = n3
            if np.linalg.norm(l.n1.X - l.normal - n.X) > np.linalg.norm(l.n1.X + l.normal - n.X):
                l.invert_normal()

            # area along x-axis: circle ring
            r_in = min([l.n1.x, l.n2.x])
            r_out = max([l.n1.x, l.n2.x])
            A_x = np.pi * (r_out**2 - r_in**2)
            # area along y-axis: cylinder
            dy = abs(l.n2.y - l.n1.y)
            r_mean = np.mean([l.n1.x, l.n2.x])
            A_y = 2 * np.pi * r_mean * dy

            center = l.n1.X + 0.5*(l.n2.X - l.n1.X)
            T_grad = t.B_e(center[0], center[1]) @ t.T

            Q += - lmbd * (A_y * T_grad[0] * np.sign(l.normal[0]) + A_x * T_grad[1] * np.sign(l.normal[1]))
            x_quiver.append(center)
            if not normal_proj:
                q_quiver.append(-lmbd * T_grad)
            else:
                q_quiver.append(-lmbd * T_grad @ l.normal * l.normal)
        fluxes.update({boundary: float(Q)})

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    x_quiver = np.array(x_quiver)
    q_quiver = np.array(q_quiver)
    ax.quiver(x_quiver[:, 0], x_quiver[:, 1], q_quiver[:, 0], q_quiver[:, 1], width=2e-3)
    ax.axis('equal')
    if save:
        with open(sim_path + '/results/heat-fluxes.yml', 'w') as f:
            yaml.dump(fluxes, f)
        plt.savefig(sim_path + '/results/heat-fluxes.pdf')
    if plot:
        plt.show()
    return fig, ax, fluxes


def boundary_scalars(sim_path):
    df_boundaries = pd.read_table(sim_path + '/boundaries.txt', delim_whitespace=True)
    boundaries = df_boundaries['Boundary-Name'].to_list()
    names_data = read_names_file(sim_path + '/results/boundary_scalars.dat.names')
    res = []
    for element in names_data:
        if 'res: ' in element:
            res.append(element[5:])
    header = []
    for flux in ['T-loads', 'Q-flux']:
        for boundary in boundaries:
            header.append(flux + '_' + boundary)
    header += ['min_flux', 'max_flux'] + res
    df = pd.read_table(sim_path + '/results/boundary_scalars.dat', names=header, sep= ' ', skipinitialspace=True)
    for column in header:
        if 'T-loads' in column or column == 'eddy current power':
            df[column] *= 2 * np.pi
    df.to_csv(sim_path + '/results/boundary_scalars.csv', index=False, sep=';')


def probes(sim_path):
    with open(sim_path + '/probes.txt', 'r') as f:
        probes = f.read().splitlines()
    names_data = read_names_file(sim_path + '/results/probes.dat.names')
    values = []
    res = []
    for column in names_data:
        if 'value: ' in column:
            values.append(column[7:].split(' in element ')[0])  # remove 'value: ' and 'in element XX'
        if 'res: ' in column:
            res.append(column[5:])
    values = values[:int(len(values)/len(probes))]  # remove duplicates
    header = []
    for probe in probes:
        for value in values:
            header.append(value + ' ' + probe)
    header_shortened = []
    for column in header:
        if 'temperature' in column and not 'loads' in column:
            header_shortened.append(column)
    header += res
    header_shortened += res
    df = pd.read_table(sim_path + '/results/probes.dat', names=header, sep=' ', skipinitialspace=True)
    df = df[header_shortened]
    if 'eddy current power' in header:
        df['eddy current power'] *= 2*np.pi
    df.to_csv(sim_path + '/results/probes.csv', index=False, sep=';')


def boundary_lines(sim_path):
    """Deprecated. Use heat_flux instead."""
    df_boundaries = pd.read_table(sim_path + '/boundaries.txt', delim_whitespace=True)
    header = ['iteration', 'BC', 'node_idx', 'x', 'y', 'z', 'T-grad_x', 'T-grad_y', 'T-loads']
    df = pd.read_table(sim_path + '/results/boundary_lines.dat', names=header, sep= ' ', skipinitialspace=True)
    df = df[df['iteration'] == max(df['iteration'])]  # take last iteration only
    # add column with bc name
    df['BC_name'] = ''
    for i, line in df_boundaries.iterrows():
        df.loc[df['BC'] == line['ID'], ['BC_name']] = line['Boundary-Name']
    boundary_conditions = df['BC_name'].unique()
    # initialize dicts for evaluation
    T_flux_dict = {bc: 0 for bc in boundary_conditions}
    T_load_dict = {bc: 0 for bc in boundary_conditions}
    # iterate over BCs
    for bc in boundary_conditions:
        df_bc = df.loc[df['BC_name'] == bc]
        # sort dict to get nodes in mathematical positive order
        # (works for rectangles or lines parallel to coordinates only)
        x_min = df_bc['x'].min()
        x_max = df_bc['x'].max()
        y_min = df_bc['y'].min()
        y_max = df_bc['y'].max()
        df_parts = []
        bottom = df_bc.loc[df_bc['y'] == y_min].sort_values(by='x')
        if len(bottom.index) > 2:
            df_parts.append(bottom)
        right = df_bc.loc[df_bc['x'] == x_max].sort_values(by='y')
        if len(right.index) > 2:
            df_parts.append(right)
        top = df_bc.loc[df_bc['y'] == y_max].sort_values(by='x', ascending=False)
        if len(top.index) > 2 and y_min != y_max:
            df_parts.append(top)
        left = df_bc.loc[df_bc['x'] == x_min].sort_values(by='y', ascending=False)
        if len(left.index) > 2 and x_min != x_max:
            df_parts.append(left)
        if df_parts != []:
            df_bc = pd.concat(df_parts)
        df_bc = df_bc.drop_duplicates()  # don't use inplace=True here to avoid warning
        df_bc.reset_index(inplace=True, drop=True)
        # T-loads (radiation, exact)
        T_load_dict[bc] = df_bc['T-loads'].sum() * 2 * np.pi
        # T-flux (radiation + conduction, approximate)
        # flux over area behind node (half distance to this node)
        for i in range(len(df_bc.index) - 1):
            T_flux_i = np.array([df_bc.iloc[i]['T-grad_x'], df_bc.iloc[i]['T-grad_y']])
            x_i = np.array([df_bc.iloc[i]['x'], df_bc.iloc[i]['y']])
            x_i1 = np.array([df_bc.iloc[i+1]['x'], df_bc.iloc[i+1]['y']])
            dx = x_i1[0] - x_i[0]
            dy = x_i1[1] - x_i[1]
            dA_x = 2 * np.pi * x_i[0] * abs(dy) / 2  # area in x-direction: cylinder
            if dx > 0:
                dA_y = np.pi * ((x_i[0] + dx/2)**2 - x_i[0]**2)  # area in y-direction: circle ring
                e_y = np.array([0, 1])
            else:
                dA_y = np.pi * (x_i[0]**2 - (x_i[0] + dx/2)**2)
                e_y = np.array([0, -1])
            if dy > 0:
                e_x = np.array([-1, 0])
            else:
                e_x = np.array([1, 0])
            T_flux_dict[bc] += np.dot(T_flux_i, e_x) * dA_x + np.dot(T_flux_i, e_y) * dA_y
        # flux over are in front of node (half distance to this node)
        for i in range(len(df_bc.index) - 1):
            T_flux_i1 = np.array([df_bc.iloc[i+1]['T-grad_x'], df_bc.iloc[i+1]['T-grad_y']])
            x_i = np.array([df_bc.iloc[i]['x'], df_bc.iloc[i]['y']])
            x_i1 = np.array([df_bc.iloc[i+1]['x'], df_bc.iloc[i+1]['y']])
            dx = x_i1[0] - x_i[0]
            dy = x_i1[1] - x_i[1]
            dA_x = 2 * np.pi * x_i1[0] * abs(dy) / 2
            if dx > 0:
                dA_y = np.pi * (x_i1[0]**2 - (x_i1[0] - dx/2)**2)
                e_y = np.array([0, 1])
            else:
                dA_y = np.pi * ((x_i1[0] - dx/2)**2 - x_i1[0]**2)
                e_y = np.array([0, -1])
            if dy > 0:
                e_x = np.array([-1, 0])
            else:
                e_x = np.array([1, 0])
            T_flux_dict[bc] += np.dot(T_flux_i1, e_x) * dA_x + np.dot(T_flux_i1, e_y) * dA_y
    # write to file
    header = ''
    values = ''
    for key, value in T_load_dict.items():
        header += 'T_loads-' + key + ';'
        values += str(value) + ';'
    for key, value in T_flux_dict.items():
        header += 'T_grad-' + key + ';'
        values += str(value) + ';'
    with open(sim_path + '/results/boundary_lines.csv', 'w') as f:
        f.write(header + '\n')
        f.write(values + '\n')


def read_names_file(names_file, skip_rows=7):
    with open(names_file, 'r') as f:
        data = f.readlines()
    data = data[7:]  # remove header
    for i in range(len(data)):
        data[i] = data[i][6:-1]  # remove index and \n
    return data


def postprocessing(simulations_file):
    with open(simulations_file, 'r') as f:
        simulations = yaml.safe_load(f)
    for simulation in simulations:
        probes(simulations[simulation])
        boundary_scalars(simulations[simulation])
        heat_flux(simulations[simulation])


if __name__ == "__main__":
    simulations_file = './simdata/simulations.yml'
    postprocessing(simulations_file)
