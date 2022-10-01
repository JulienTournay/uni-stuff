import numpy as num
import sympy as sym
import statistics as stats

tab = '    '
begin = '{'
end = '}'
number = '\\num{'


def standardabweichung(list): return num.sqrt(num.sum([(i-stats.mean(list))**2 for i in list])/(len(list)-1))


def get_data_csv(csv_path): return [[line[i] for line in num.recfromcsv(csv_path, delimiter=',')] for i in range(len(num.recfromcsv(csv_path, delimiter=',')[1]))]


def get_mean_csv(csv_path): return [stats.mean(i) for i in [[line[i] for line in num.recfromcsv(csv_path, delimiter=',')] for i in range(len(num.recfromcsv(csv_path, delimiter=',')[1]))]]


def get_mean(data): return [stats.mean(i) for i in data]


# PRINT A USABLE LATEX TABLE DIRECTLY FROM A .CSV FILE.
# csv_path: Path to the.csv file that is to be made into a table.
# columns: List of ints of the columns to be printed starting at 0. If left blank, all columns will be printed.
# rows: List of ints of the rows to be printed starting at 0. If left blank, all rows will be printed.
def make_table_csv(csv_path, columns=False, rows=False):
    data = get_data_csv(csv_path)
    name = csv_path.split('/')[-1].split('.')[0]

    if not columns: columns = list(range(len(data)))
    clmns = [data[n] for n in columns]

    if not rows: rows = list(range(len(data[0])))

    head = f'\n{tab}'.join([
        '\\begin{table}[]',
        '\\centering',
        '\\begin{tabular}{|' + 'c|'*len(columns) + '}',
        tab + '\\hline'
        ])

    foot = tab*2 + f'\n{tab}'.join([
        '\\hline',
        '\\end{tabular}',
        '\\caption{' + name + '}',
        '\\label{tab:' + name + '}\n\\end{table}'
        ])

    print(head)
    for n in rows:
        cont = f'{tab*2}'
        first = True
        for item in clmns:
            if first:
                cont += f'{number}{item[n]}{end}'
                first = False
            else: cont += f' & {number}{item[n]}{end}'
        cont += ' \\\\'
        print(cont)
    print(foot)


# To be enhanced.
def make_table(data):
    print('\\begin{table}[]\n    \\centering\n    \\begin{tabular}{|' + ''.join(['c|' for item in data]) + '}\n        \\hline')
    for n in range(len(data[0])):
        cont = f'{tab*2}'
        first = True
        for item in data:
            if first:
                if isinstance(item[n], num.number): cont += f'{number}{item[n]:3.2f}{end}'
                elif isinstance(item[n], sym.Float): cont += f'{number}{item[n]:3.2f}{end}'
                else: cont += f'{item[n]}'
                first = False
            else:
                if isinstance(item[n], num.number): cont += f' & {number}{item[n]:3.2f}{end}'
                elif isinstance(item[n], sym.Float): cont += f' & {number}{item[n]:3.2f}{end}'
                else: cont += f' & {item[n]}'
        cont += ' \\\\'
        print(cont)
    print("        \\hline\n    \\end{tabular}\n    \\caption{Caption}\n    \\label{tab:my_label}\n\\end{table}")


def gauss(formula, vars, space=1, other_vars=[]):
    new_stack = []
    for n in range(space):
        list_formeln = []
        for i in range(len(vars)):
            derivative = formula.diff(vars[i][0])
            for j in vars: derivative = derivative.subs(j[0], j[1][n])
            for k in other_vars: derivative = derivative.subs(k[0], k[1])
            list_formeln.append((derivative * vars[i][2])**2)
        new_stack.append(sym.Float(sym.sqrt(sum(list_formeln))))
    return new_stack


def gauss_latex(formula, vars):
    list_formulas = []
    for i in range(len(vars)):
        derivative = formula.diff(vars[i][0])
        list_formulas.append(f'&+ \\left({sym.latex(derivative)} \\Delta {sym.latex(vars[i][0])} \\right)^2 \\nonumber \\\\')
    print('\\begin{align}')
    for n in range(len(list_formulas)): print(tab + list_formulas[n])
    print(tab + '\\label{gaussformel}\n\\end{align}')