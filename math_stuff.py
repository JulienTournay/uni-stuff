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


def get_mean_list(list): return [stats.mean(i) for i in list]


# Rotate a list by 90 degrees. List must be square.
def rotate_list(list): return [[item[i] for item in list] for i in range(len(list[0]))]


# PRINT A USABLE LATEX TABLE DIRECTLY FROM A .CSV FILE.
# csv_path: Path to the .csv file that is to be made into a table.
# columns: List of ints of the columns to be printed starting at 0. If left blank, all columns will be printed.
# rows: List of ints of the rows to be printed starting at 0. If left blank, all rows will be printed.
# rotate: Bool. Will rotate table by 90 degrees if True.
# from: List of length 2. First entry: format (e.g. .2 will put out '100.00' for an input of 100 or '23.98' for an input of 23.98658). Second entry: format type (e.g. 'f').
def make_table_csv(csv_path, columns=False, rows=False, rotate=False, form=False):
    data = get_data_csv(csv_path)
    name = csv_path.split('/')[-1].split('.')[0]

    if rotate: data = rotate_list(data)

    if not columns: columns = list(range(len(data)))
    clmns = [data[n] for n in columns]

    if not rows: rows = list(range(len(data[0])))

    if not form: form = [.2, 'f']

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
                cont += f'{number}{"{:{}{}}".format(item[n], form[0], form[1])}{end}'
                first = False
            else: cont += f' & {number}{"{:{}{}}".format(item[n], form[0], form[1])}{end}'
        cont += ' \\\\'
        print(cont)
    print(foot)


# PRINT A USABLE LATEX TABLE FROM A LIST.
# csv_path: Path to the .csv file that is to be made into a table.
# columns: List of ints of the columns to be printed starting at 0. If left blank, all columns will be printed.
# rows: List of ints of the rows to be printed starting at 0. If left blank, all rows will be printed.
# rotate: Bool. Will rotate table by 90 degrees if True.
# from: List of length 2. First entry: format (e.g. .2 will put out '100.00' for an input of 100 or '23.98' for an input of 23.98658). Second entry: format type (e.g. 'f').
def make_table_list(data, columns=False, rows=False, rotate=False, form=False):

    if rotate: data = rotate_list(data)

    if not columns: columns = list(range(len(data)))
    clmns = [data[n] for n in columns]

    if not rows: rows = list(range(len(data[0])))

    if not form: form = [3.2, 'f']

    def instancecheck(var):
        if isinstance(var, str): return var
        else: return f'{number}{"{:{}{}}".format(var, form[0], form[1])}{end}'

    head = f'\n{tab}'.join([
        '\\begin{table}[]',
        '\\centering',
        '\\begin{tabular}{|' + 'c|'*len(columns) + '}',
        tab + '\\hline'
        ])

    foot = tab*2 + f'\n{tab}'.join([
        '\\hline',
        '\\end{tabular}',
        '\\caption{name}',
        '\\label{tab:name}\n\\end{table}'
        ])

    print(head)
    for n in rows:
        cont = f'{tab*2}'
        first = True
        for item in clmns:
            if first:
                cont += instancecheck(item[n])
                first = False
            else: cont += f' & {instancecheck(item[n])}'
        cont += ' \\\\'
        print(cont)
    print(foot)


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