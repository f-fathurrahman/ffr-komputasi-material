# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3.9.13 ('base')
#     language: python
#     name: python3
# ---

# %%
import numpy as np
from sympy import symbols, Symbol, diff, log
import re

output_style = 'taichi' # numpy or taichi
supplementary_infomation = False # Print supplementary information?｜補足情報をprintするか
useful_function_for_PFM = True # Print function for phase-field simulation?｜PF計算用の関数をprintするか
userful_function_for_PFM_Ni_superalloy = False
full_TDB_text = False # Print original TDB file at the end?｜最後に元のTDBファイルを出力するか
R = 8.3145 # Gus constant｜気体定数
Tfix_True = False # これは計算速度にあまり影響しないから、いらない
Tfix = 1400.0 # 等温保持計算用



# file_path = '../tdb_file/mc_fe_fcc.tdb'
# solvent   = 'FE'
# # solutes   = ['CR', 'NI'] #溶質
# # solutes   = ['AL', 'CO', 'CR', 'CU', 'NI']
# # solutes   = ['AL', 'CO', 'CR', 'CU', 'HF', 'LA', 'MN', 'NI']
# # solutes   = ['AL', 'CO', 'CR', 'CU', 'HF', 'LA', 'MN', 'MO', 'NB', 'NI', 'P']
# solutes   = ['AL', 'CO', 'CR', 'CU', 'HF', 'LA', 'MN', 'MO', 'NB', 'NI', 'P', 'PD', 'S', 'SI']
# # solutes   = ['AL', 'CO', 'CR', 'CU', 'HF', 'LA', 'MN', 'MO', 'NB', 'NI', 'P', 'PD', 'S', 'SI', 'TA', 'TI', 'V']
# # solutes   = ['AL', 'CO', 'CR', 'CU', 'HF', 'LA', 'MN', 'MO', 'NB', 'NI', 'P', 'PD', 'S', 'SI', 'TA', 'TI', 'V', 'W', 'Y'] #溶質
# phases = ['LIQUID', 'FCC_A1']


file_path = '../tdb_file/Ni_superalloy.tdb'
phases = ['FCC_A1', 'L12']
solvent   = 'NI'
solutes = ['AL', 'CO'] #溶質
# solutes = ['AL', 'CO', 'CR', 'CU', 'FE'] #溶質
# solutes = ['AL', 'CO', 'CR', 'CU', 'FE', 'HF', 'MN', 'MO'] #溶質
# solutes = ['AL', 'CO', 'CR', 'CU', 'FE', 'HF', 'MN', 'MO', 'NB', 'SI', 'TI'] #溶質
# # solutes = ['AL', 'CO', 'CR', 'CU', 'FE', 'HF', 'MN', 'MO', 'NB', 'SI', 'TI', 'V', 'W', 'ZR'] #溶質

# file_path = '../tdb_file/Ni_superalloy.tdb'
# solvent   = 'NI'
# # solutes   = ['AL', 'CO'] #溶質
# solutes     = ['AL', 'CO', 'CR', 'FE', 'HF'] #溶質
# # solutes     = ['AL', 'CO', 'CR', 'FE', 'HF', 'MO', 'NB', 'SI'] #溶質
# # solutes     = ['AL', 'CO', 'CR', 'FE', 'HF', 'MO', 'NB', 'SI', 'TI', 'W', 'ZR'] #溶質
# phases = ['FCC_A1', 'L12']

# file_path = '../tdb_file/mc_al_fcc.tdb'
# phases = ['LIQUID', 'FCC_A1']
# solvent = 'AL'
# # solutes = ['CU','MG']
# # solutes = ['CR', 'CU', 'FE', 'MG', 'MN']
# solutes = ['CR', 'CU', 'FE', 'MG', 'MN', 'NI', 'SC', 'SI']
# # solutes = ['CR', 'CU', 'FE', 'MG', 'MN', 'NI', 'SC', 'SI', 'TI', 'ZN', 'ZR']
# # solutes = ['CR', 'CU', 'FE', 'MG', 'MN', 'SI', 'TI', 'ZN', 'ZR'] # AA7050

# %% [markdown]
# # Capital letter variables
# - LINES （list devided by "!"）
# - LINES_FUNCTION (list of FUNCTION)
# - LINES_CONSTITUENT (list of CONSTITUENT)
# - LINES_PARAMETER (list of PARAMETER)
# - LINES_PHASE (list of PHASE)
# - CONFIGURATION (Dictionary of each phase→ ['site'] (stoichiometry of each sublattice) ['element'] (element of each sublattice))
# - FUNCTIONS ['FUNCTION name'] (text of each FUNCTION)
#

# %%
with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
    LINES = f.readlines() # Read as list (with EOL)｜リストとして読み込む（改行コードあり）


# Adjust LINES｜LINES調整ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
LINES = [line.strip('\n') for line in LINES] # Remove EOL｜改行コードを取り除く
LINES = [line.strip() for line in LINES] # Remove spaces in the beginning and end of sentence｜ 文頭文末から空白を取り除く
LINES = [line for line in LINES if not line.startswith('$')] # Remove comment out lines｜コメントアウト行を削除
LINES = [line+' ' for line in LINES] # Insert space into the beginning of sentence｜文頭にスペース挿入
LINES = ''.join(LINES) # Join all lists｜リストを全て結合
LINES = LINES.replace('#', '') # Remove #｜#を削除する
LINES = LINES.upper() # Convert to capitals｜大文字にする

# Remove option embedded with phase name｜相の名前に付くオプションの削除
options = ['L','G','A','Y','I','F','B','C'] # L:liquid, G:gus, A:aqueous solution, ... etc.
for option in options:
    for input_phase in phases:
        pattern = input_phase + ':' + option
        LINES = LINES.replace(pattern, input_phase) # Remove option from (input_phase+option)｜input_phase+オプションの文字列からオプション部分を削除する

# Adjust LINES｜LINES調整
LINES = LINES.split('!') # Make list by separationg "!"｜!で区切ってリストを作成
LINES = [re.split('>', line)[0] for line in LINES] # Remove after ">"｜> 以降を削除する
LINES = [re.split(' REF:', line)[0] for line in LINES] # Remove after "REF:"｜REF:以降を削除する
LINES = [line.strip() for line in LINES] # Remove spaces in the beginning and end of sentence｜ 文頭文末から空白を取り除く

# Make capital letter variables
LINES_FUNCTION    = [line for line in LINES if line.startswith('FU')] # LINES of FUNCTION｜FUNCTION部分
LINES_CONSTITUENT = [line for line in LINES if line.startswith('CO')] # LINES of CONSTITUENT｜CONSTITUENT部分
LINES_PARAMTER    = [line for line in LINES if line.startswith('PA')] # LINES of PARAMTER｜PARAMTER部分
LINES_PHASE       = [line for line in LINES if line.startswith('PH')] # LINES of PHASE｜PHASE部分



# Make CONFIGURATION｜CONFIGURATION作成ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
CONFIGURATION = {}
for input_phase in phases:

    # Retrieve site information｜site情報の抽出
    PH = [item for item in LINES_PHASE if input_phase in item][0] # Retrieve PHASE line including input_phase｜input_phaseのPHASEラインを取得
    pattern = r'\s+\d\s+' # Pattern of "space+number+space"｜スペース+数字一文字+スペース
    match = re.search(pattern, PH)
    PH = PH[match.start():] # Retrieve the searched part and onwards｜検索した部分とそれ以降を取得
    PH = PH.split() # Devide by space｜空白で区切る
    number_of_sublattice = int(PH[0]) # number_of_sublattice｜副格子の数
    site_list = [] 
    for sublattice in range(number_of_sublattice):
        site_list.append(float(PH[1+sublattice])) # Retrieve stoichiometry of each sublattice｜各副格子のサイト占有率をsite_listに入れる

    # Retrieve element information｜element情報の抽出
    CO = [item for item in LINES_CONSTITUENT if input_phase in item][0] # Retrive CONSTITUENT line including input_phase｜input_phaseが入っているCONSTITUENTラインを取得
    CO = CO.replace("%","").replace(" ","") # Remove "%" and space｜%とスペースを削除
    CO = CO.split(':')[1:-1] # Make list devided by ":" and ignore first and last component｜:で区切ってリスト作成し、リストの最初と最後の成分を無視
    CO = [part.split(',') for part in CO] # Devide each list by ","｜リストの各成分のテキストを,で区切る
    CO = [[element for element in sublist if element in [solvent]+solutes] for sublist in CO] # Retrieve element that correspond to solvent or solute｜ solvent or soluteに該当するelementのみ抽出

    # Remove sublattices that do not need to be calculated｜今回計算する必要のない副格子を削除する
    element_list = []
    for i,ele in enumerate(CO):
        if not ele == []: # 空のリストでない場合はelement_listに追加する
            element_list.append(ele)
        else: # 空のリストの場合、その副格子を削除する
            del site_list[i] 

    CONFIGURATION[input_phase] = {'site':site_list, 'element':element_list}



# Make FUNCTIONS｜FUNCTIONS作成ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
FUNCTIONS = {}
for FU in LINES_FUNCTION.copy():
    FUNCTION_name = (FU.split())[1] # Name of FUCTION｜FUNCTIONの関数名
    FUNCTIONS[FUNCTION_name] = None # Add key｜キーを追加

# Function to extract the function name and function from the FUNCTION (or PARAMETER) line｜FUNCTION (or PARAMETER) ラインから関数名と関数を抜き出す関数
def extract_function_name_function(input_line, PARAMTER_line=False): 
    T_domain = [] # List of temperature range｜温度域のリスト
    formula = [] # List of formula for each temperature range｜各温度域での式のリスト
    FU = input_line.split(' Y') # Devide each line by "Y"｜各行を Yで区切る

    for i,FU_Y in enumerate(FU): # Make function for each temperature range｜各温度帯の関数を作成する
        pattern = r'\sN\s|\sN!|\sN$' # Pattern of "space+N+space" or "space+N+!" or "space+N+$" ｜パターン（スペース + N + スペース or スペース + N + ! or スペース+N+文末($)）
        FU_Y = re.split(pattern, FU_Y)[0] # Delete after the pattern｜パターンにマッチした部分以降を削除
        FU_Y = FU_Y.split() #　Devide by space｜空白で区切る
        if i == 0: # First temperature range｜最初の温度定義
            func_name = FU_Y[1] # Name of function｜関数の名前

            if PARAMTER_line: # Give my own function names when dealing with PARAMTER lines with temperature branches｜温度分岐のあるPARAMTERラインを取り扱う場合、独自に関数名をつける
                pattern = r'[(),;: ]'
                func_name = re.sub(pattern, '', func_name) # Replace "(),;:" to "_" ｜ (),;:スペースを全て_に置換する

            T_domain.append(FU_Y[2]) # Initial defined temperature｜定義開始温度
            FU_Y = FU_Y[3:] # Delete "FUNCTION", function name and initial defined temperature｜ "FUNCTION", "関数名", "定義開始温度" を削除

        T_domain.append(FU_Y[-1]) # Last defined temperature｜定義終了温度
        FU_Y.pop(-1) # Delete last defined temperature｜定義終了温度を削除
        FU_Y = ''.join(FU_Y) # Join all list to address newline｜リストを結合（改行に対応するため）
        FU_Y = FU_Y.replace(';','').replace('LN(T)','ti.log(T)').replace('ln(T)','ti.log(T)').replace('T**(-','(1/T)**(')
        for NAME in list(FUNCTIONS.keys()): # Add (T) to FUNCTION｜FUCNTIONに(T)をつける
            pattern = rf'(?<=[^\w])({NAME})(?=[^\w])|(?<=[^\w])({NAME})$' # Pattern of "not alphabet + FUNCTION name + not alphabet" or "not alphabet + FUNCTION name + $"｜アルファベット以外+FUNCTION name+アルファベット以外　or アルファベット以外+FUNCTION name+文末
            FU_Y = re.sub(pattern, r'\g<0>(T)', FU_Y) # Add (T)｜(T)を加える
        formula.append(FU_Y)

    # Construct fuction code｜functionのコードを作成する
    function = f'@ti.func\ndef {func_name}(T):\n    x = 0.0\n    if T >= {T_domain[0]} and T < {T_domain[1]}:\n        x =  {formula[0]}\n'
    for i in range(1,len(formula)):
        function += f'    elif T >= {T_domain[i]} and T < {T_domain[i+1]}:\n        x =  {formula[i]}\n'
    function += '    return x\n'

    return func_name, function


# Make FUNCTIONS｜FUNCTIONSの作成
for FU in LINES_FUNCTION:
    func_name, function = extract_function_name_function(FU, PARAMTER_line=False)
    FUNCTIONS[func_name] = function


# %% [markdown]
# # formulate
# - Make free energy functions from PARAMTER lines｜PARAMETERラインから、自由エネルギー式を作成する

# %%
def formulate(input_phase, info:False):

    def PARAMETER_infomation(PARAMETER_text): # Extract function, index n, 成分｜PARAMETERライン（一行）から、関数、指数ｎ、成分を抜き出す
        function = PARAMETER_text
        pattern = r";\s*\d\s*\)(.*)" # "; + space + number + space + )"｜; + 'スペースok' + '数字' + 'スペースok' + ) 
        match = re.search(pattern, function)
        function = match.group(1) # Get text after ")"｜)より後を取得

        function = function.strip() # Remove spaces at the beginning and end of text｜頭尾のスペースを消す
        function = function.split(' ') # Devide by space｜スペースで切る
        function = function[1:] # Delete initial defined temperature｜定義開始温度部分を削除
        function = ''.join(function)
        function = function.split(';') # Devide by ;｜;で区切る
        function = function[:-1] # ; Delete end defined temperature｜定義終了温度を削除
        function = ''.join(function)
        function = '(' + function + ')'

        pattern = r";(\d)\)" # "; + number + )"　｜　;+数字+)
        match = re.search(pattern, PARAMETER_text) 
        index = match.group(1).strip() # index n｜ インデックス部分(指数ｎのこと)

        pattern = pattern = r"[A-Z]\(([^,]+),([^;]+)\;" #  pattern of G(...)｜G(...)部分
        match = re.search(pattern, PARAMETER_text) 
        component = match.group(2).strip()
        component = re.split('[:]',component) # Devide by : ｜ :で区切って副格子ごとにする
        for i,comp in enumerate(component): # Devide by , in each sublaattice｜各副格子の中の , で区切る
            component[i] = re.split('[,]',comp)
        return function, index, component
    
    def check_membership(A, B): # Determine whether all components of B (2D ndarray) are included in A (1D list)｜B(2次元ndarray)の全ての成分がA(１次元リスト)に含まれているかの判定
        A_set = set(A)
        for row in B:
            for element in row:
                if element not in A_set:
                    return False
        return True
    
    def diploy_wild_card(PARAMETERS): # Wildcard expansion (incomplete→Not available if there are two * in one line.)｜ワイルドカードの展開（未完成→一つの行に二つの＊が入っている場合は無理）
        index_slide = 0 # Correction of index deviations by substitution｜置換によるindexのずれの補正
        for i, PARAMETER in enumerate(PARAMETERS.copy()):
            function, index, component = PARAMETER_infomation(PARAMETER)
            for j, component_sublattice in enumerate(component): # component for each sublattice｜各副格子のcomponent
                if '*' in component_sublattice: # Exist * in the sublattice?｜この副格子にワイルドカードが入っているかどうか判定
                    if info:
                        print(f'# diployed wild card in {component}')
                    replace_comps = CONFIGURATION[input_phase]['element'][j] # Components to which wildcards are substituted｜ワイルドカードを置換させる先の成分
                    replace = [] 
                    for replace_comp in replace_comps:
                        replaced_text = PARAMETER.replace('*', replace_comp, 1) # Replace the * that appears first｜一番初めに登場する*のみを置換する→ここが未完成→一つの行に二つの＊が入っている場合は無理
                        replace.append(replaced_text)
                    slided_i = i + index_slide
                    del PARAMETERS[slided_i] # Delete line including *｜ワイルドカードの行を削除
                    PARAMETERS[slided_i:slided_i] = replace # Insert the expanded list at the deleted location｜削除した場所にワイルドカード展開後のリストを挿入する
                    index_slide += len(replace)-1 # Adjust index　displacement due to the insert｜挿入するとindexがずれるので調整する
        return PARAMETERS
    
    def formulate_entropy(): # Make entropy function｜エントロピー項の式を作成する
        elements = CONFIGURATION[input_phase]['element']
        sites    = CONFIGURATION[input_phase]['site']
        site_sum = sum(sites)
        phase_formula = f'1/{site_sum}*(+ \\\n           + {R}*T*('
        for i, elements in enumerate(elements):
            phase_formula += f' +{sites[i]}*('
            for X in elements:
                X = X+str(i+1) # Add sublattice number to element｜elementに副格子番号つける
                phase_formula += f'+({X})*ti.log({X}) ' 
            phase_formula += ') '
        phase_formula += ')'
        return phase_formula
    
    def check_branch(PARAMTER): # Check if there is temperature branch ｜温度分岐があるかどうかをチェック
        PARAMTER = PARAMTER.split(';') # Devide by ;｜;で分割
        if len(PARAMTER) > 3:
            return True
        else:
            return False
        

    phase_formula = formulate_entropy() # Make function of entropy｜エントロピー項のエネルギー式を作成する
    number_of_sublattice = len(CONFIGURATION[input_phase]['site'])
    PARAMETERS = [item for item in LINES_PARAMTER if input_phase in item] # Get PHASE line including input_phase｜input_phaseのPHASEラインを取得
    PARAMETERS = [item for item in PARAMETERS if re.search(r' G\(| L\(', item)] # Get line including G( or L(｜" G(" or " L("がある行を抽出
    PARAMETERS = diploy_wild_card(PARAMETERS) # Diploy wild cards｜ワイルドカードを展開する

    for i, PARAMETER in enumerate(PARAMETERS):
        function, index, component = PARAMETER_infomation(PARAMETER)

        # Define FUNCTION for temperature branch｜温度分岐部分をFUNCTION定義する
        temperature_branch = check_branch(PARAMETER)
        if temperature_branch:
            func_name, function = extract_function_name_function(PARAMETER, PARAMTER_line=True)
            FUNCTIONS[func_name] = function # Add function to the dictionary of FUCNTION｜FUNCTIONの辞書に追加
            function = func_name # Replace function by func_name｜functionをfunc_nameに置き換える必要がある


        valid_endmember = check_membership([solvent]+solutes+['VA'], component) 
        if valid_endmember: # Check if the endmember consists of the component you want to calculate｜計算したい成分で構成されたエンドメンバーかを確認

            if info:
                print(f'# {component} {index} : {function}')
        
            if component[-1][0] == 'VA': # If the last sublattice is 'VA' only, delete it｜ 最後の副格子が'VA'単体なら削除する
                component = [sublattice for sublattice in component if sublattice != ['VA']]
                number_of_sublattice += -1 # Decrease the number of sublattice｜副格子数を減らす

            coefficient = '' # Composition part of the thermodynamic function｜ 関数の濃度部分
            for i,comp in enumerate(component):
                comp =  [X+f'{i+1}' for X in comp] # Add sublattice number｜副格子番号を付ける
                comp =  [X+' ' if len(X) ==2 else X for X in comp] # Insert space to make it easier to read｜適宜スペースを挿入して見やすくする
                    
                if   len(comp) == 1: # Urnary｜純物質
                    A = f'{comp[0]}                                       ' 
                
                elif len(comp) == 2: # Binary｜２元相互作用パラメータ
                    if index == '0':   A = f'{comp[0]}*{comp[1]}                                   '     
                    elif index == '1': A = f'{comp[0]}*{comp[1]}*({comp[0]}-{comp[1]})                         '     
                    else:              A = f'{comp[0]}*{comp[1]}*({comp[0]}-{comp[1]})**{index}                      '     
                
                elif len(comp) == 3: # Ternary（This code is a little bit complex to accommodate abbreviations)｜３元相互作用パラメータ（省略表記に対応させるために少し複雑になってる）
                    if int(index) == 0: # Check if this line is abbreviation or not（Search for lines consist of the same sublattice configuration as present line)｜省略表記かどうか判定していく（同じ副格子構成のラインがあるか探す）
                        omit = True
                        for PARAMETER_2 in PARAMETERS: # Check all PARAMTER lines｜全てのPARAMETERライン回る
                            function_2, index_2, component_2 = PARAMETER_infomation(PARAMETER_2)
                            if component_2[-1][0] == 'VA': # If the last sublattice is 'VA' only, delete it｜ 最後の副格子が'VA'単体なら削除する
                                component_2 = [sublattice for sublattice in component_2 if sublattice != ['VA']]
                            if component == component_2 and int(index_2) > 1: # If the configuration of sublattice is same, we don't use abbreviation｜同じ元素構成でn=1 or n=2が存在する場合は省略表記しない
                                omit = False
                        if info:
                            print(f'# ternary omit={omit}')
                        if omit:
                            A = f'{comp[0]}*{comp[1]}*{comp[2]}                               '  # Abbreviation｜省略表記
                        else:
                            A = f'{comp[0]}*{comp[1]}*{comp[2]}*({comp[0]}+(1/3)*(1-({comp[0]}+{comp[1]}+{comp[2]}))) ' # Normal notation｜通常表記
                    else: # Abbreviations are not allowed for n=1 and n=2｜n=1とn=2の場合は省略表記不可
                        A = f'{comp[0]}*{comp[1]}*{comp[2]}*({comp[int(index)]}+(1/3)*(1-({comp[0]}+{comp[1]}+{comp[2]}))) ' # Normal notation｜通常表記
                
                elif len(comp) == 4: # Quaternary｜４元相互作用パラメータ
                    A = f'{comp[0]}*{comp[1]}*{comp[2]}*{comp[3]} '
                
                coefficient += A + '*'

            phase_formula += '+ \\\n           + ' + coefficient + function # Backslash can be escaped by doubling it｜バックスラッシュは２重にすることでエスケープできる

    # Final adjustment｜最終調整
    for FUNCTION_name in list(FUNCTIONS.keys()):
        pattern = rf'(?<=[^\w])({FUNCTION_name})(?=[^\w])' # "not alphabet + FUCNTION_name + not alphabet"｜アルファベット以外+FUNCTION_name+アルファベット以外
        phase_formula = re.sub(pattern, r'\g<0>(T)', phase_formula) # Add argment (T)｜引数(T)を持たせる
        phase_formula = phase_formula.replace('LN','ti.log').replace('ln','ti.log').replace('T**(-','(1/T)**(')
    phase_formula += ' )'

    return phase_formula


# %% [markdown]
# # construct_f_fy_fyy
# - Algebraic differentiation of free energy expressions｜自由エネルギー式の代数微分

# %%
def define_FUNCTION_as_variable(Tfix):
    convert_FUNCTION_into_value = {}
    for FUN_name in FUNCTIONS.keys():
        define_FUNCTION = FUNCTIONS[FUN_name].replace('@ti.func', '').replace('ti.', 'np.')
        exec(define_FUNCTION, globals()) # FUNCTIONをglobalでdefineする
        fixed_value = eval(FUN_name+f'({Tfix})') # FUNCTIONにTfixを入れた時の値
        # def_variable = FUN_name + ' = ' + str(fixed_value) # 変数として宣言するtext
        # def_variables.append(def_variable)
        convert_FUNCTION_into_value[FUN_name + '(T)'] = str(fixed_value)

    # for def_variable in def_variables:
    #     exec(def_variable, globals())
    #     # print(def_variable)
    return convert_FUNCTION_into_value # dict['FUNCTINON's name] = value when Tfix



def construct_f_fy_fyy(phase_formula, input_phase, info=False):

    
    phase_formula_sympy = phase_formula # text for differentiation by sympy｜sympyで微分する用のtext

    # Rewrite solvent by solutes in each sublattice｜溶媒濃度を溶質濃度で書き換える
    text = ''
    elements = CONFIGURATION[input_phase]['element']
    for i, element in enumerate(elements):
        equal_main_component = '1.0'
        for comp in element:
            if not comp == solvent: # If solute｜溶質なら
                # X = comp + str(i+1) # Add sublattice number｜副格子番号を付ける
                equal_main_component += f'-{comp + str(i+1)}'
        equal_main_component = '('+equal_main_component+')'
        text += (f'{solvent}{str(i+1)} = {equal_main_component}\n    ') # Code to rewrite solvent by solutes｜溶媒濃度を溶質濃度に置き換えるコード
        search_text = f'([^A-Za-z]+){solvent}{str(i+1)}' # "not alphabet + solvent"｜ アルファベット以外＋'溶媒'
        substitute_text = r'\1' + equal_main_component 
        phase_formula_sympy = re.sub(search_text, substitute_text, phase_formula_sympy)

    # Preprocessing for algebraic calculations with Sympy｜Sympyで代数計算するための前処理
    phase_formula_sympy = phase_formula_sympy.replace('\\','').replace('ti.log','log')
    
    # Tfix_True = True
    if Tfix_True:
        convert_FUNCTION_into_value = define_FUNCTION_as_variable(Tfix)
        for key, value in convert_FUNCTION_into_value.items():
            phase_formula_sympy = phase_formula_sympy.replace(key, value)

        # phase_formula_sympyの式を整理する
        phase_formula_sympy = sp.sympify(phase_formula_sympy) # textのsympy化
        T = Symbol('T')
        phase_formula_sympy = phase_formula_sympy.subs({T: Tfix}) # Tの代入
    

    declared_symbols = []
    for i, element in enumerate(elements):
        for comp in element:
            if not comp == solvent: # 溶質なら
                X = str(comp + str(i+1))
                declared_symbols.append(X)
                X = Symbol(X)

    # Make argument for Taichi function｜Taichi関数の引数の作成
    argu = 'T:ti.f64'
    for X in declared_symbols:
        argu += f',{str(X)}:ti.f64'


    # Make f_txt, fy_txt, fyy_txt｜f_txt fy_txt fyy_txtの作成
    f_txt = f'@ti.func\ndef f_{input_phase}({argu}) -> ti.f64:\n    {text}return {phase_formula}\n'
    fy_txt, fyy_txt = '', ''
    if Tfix_True:
        f_txt = f'@ti.func\ndef f_{input_phase}({argu}) -> ti.f64:\n    {text}return {phase_formula_sympy}\n'
    for X in declared_symbols:
        fy_X  = str(diff(phase_formula_sympy, X)) # Differenciate by X｜Xで微分
        fyy_X = str(diff(fy_X, X)) # Differenciate by X｜さらにXで微分

        fy_X  = fy_X.replace('log','ti.log')
        fyy_X = fyy_X.replace('log','ti.log')

        fy_txt    += f'@ti.func\ndef fy_{input_phase}_{X}({argu}) -> ti.f64:\n    fy_{X} = {fy_X}\n    return fy_{X}\n'
        fyy_txt   += f'@ti.func\ndef fyy_{input_phase}_{X}({argu}) -> ti.f64:\n    fyy_{X} = {fyy_X}\n    return fyy_{X}\n'
    if Tfix_True:
        f_txt = f_txt.replace('log', 'ti.log')
    return f_txt, fy_txt, fyy_txt


# phase_formula = formulate('FCC_A1', info=False)
# f_txt, fy_txt, fyy_txt = construct_f_fy_fyy(phase_formula, 'FCC_A1', info=False)
# print(f_txt)

# %% [markdown]
# # construct_function_for_PFM (additional)

# %%
def construct_function_for_PFM(): # PFMで使いやすい関数を作成する

    cal_f_txt   = f'@ti.func\ndef cal_f(i:int,T:ti.f64,y:ti.f64) -> ti.f64:\n    f = 0.0\n'
    cal_fy_txt  = f'@ti.func\ndef cal_fy(a:int,j:int,T:ti.f64,y:ti.f64) -> ti.f64:\n    fy = 0.0\n'
    cal_fyy_txt = f'@ti.func\ndef cal_fyy(a:int,j:int,T:ti.f64,y:ti.f64) -> ti.f64:\n    fyy = 0.0\n'

    a_accumulation = 0 # 累積副格子番号

    argu_phase = {} # 各相の引数を揃える
    for input_phase in phases:
        txt = 'T'
        number_of_sublattice = len(CONFIGURATION[input_phase]['site'])
        for a in range(number_of_sublattice):
            for j in range(len(solutes)):
                txt += f',y[{a_accumulation},{j}]'
            a_accumulation += 1

        argu_phase[input_phase] = txt


    a_accumulation = 0 # 累積副格子番号
    for i,input_phase in enumerate(phases):


        # cal_f_txtの作成
        if i == 0:
            cal_f_txt += f'    if   i == {i}: '
        else:
            cal_f_txt += f'    elif i == {i}: '
        cal_f_txt += f'f = f_{input_phase}({argu_phase[input_phase]})\n'


        # input_phaseの各副格子で場合分けする
        number_of_sublattice = len(CONFIGURATION[input_phase]['site'])
        for a in range(number_of_sublattice):
            # cal_fy_txtとcal_fyy_txtの作成
            if a_accumulation == 0:
                cal_fy_txt  += f'    if   a == {a_accumulation}: # if sublattice {a+1} in {input_phase}\n'
                cal_fyy_txt += f'    if   a == {a_accumulation}: # if sublattice {a+1} in {input_phase}\n'
            else:
                cal_fy_txt  += f'    elif a == {a_accumulation}: # if sublattice {a+1} in {input_phase}\n'
                cal_fyy_txt += f'    elif a == {a_accumulation}: # if sublattice {a+1} in {input_phase}\n'
            for j in range(len(solutes)):
                if j == 0:
                    cal_fy_txt  += f'        if   j == {j}: fy = fy_{input_phase}_{solutes[j]}{a+1}({argu_phase[input_phase]})\n'
                    cal_fyy_txt += f'        if   j == {j}: fyy = fyy_{input_phase}_{solutes[j]}{a+1}({argu_phase[input_phase]})\n'
                else:
                    cal_fy_txt  += f'        elif j == {j}: fy = fy_{input_phase}_{solutes[j]}{a+1}({argu_phase[input_phase]})\n'
                    cal_fyy_txt += f'        elif j == {j}: fyy = fyy_{input_phase}_{solutes[j]}{a+1}({argu_phase[input_phase]})\n'

            a_accumulation += 1
    
    cal_f_txt   += f'    return f\n'
    cal_fy_txt  += f'    return fy\n'
    cal_fyy_txt += f'    return fyy\n'

    return cal_f_txt, cal_fy_txt, cal_fyy_txt


def construct_function_for_PFM_superalloy(): # PFMで使いやすい関数を作成する for ni super alloy

    cal_f_txt   = f'@ti.func\ndef cal_f(i:int,T:ti.f64,y:ti.f64) -> ti.f64:\n    f = 0.0\n'
    cal_fy_txt  = f'@ti.func\ndef cal_fy(a:int,j:int,T:ti.f64,y:ti.f64) -> ti.f64:\n    fy = 0.0\n'
    cal_fyy_txt = f'@ti.func\ndef cal_fyy(a:int,j:int,T:ti.f64,y:ti.f64) -> ti.f64:\n    fyy = 0.0\n'

    argu_phase = {} # 各相の引数を揃える
    txt_FCC_A1 = 'T'
    txt_L12_1  = 'T'
    txt_L12_2  = 'T'
    for j in range(len(solutes)):
        txt_FCC_A1 += f',y[0,{j}]'
        txt_L12_1  += f',y[a,{j}]'
        txt_L12_2  += f',y[a-1,{j}]'
    for j in range(len(solutes)):
        txt_L12_1  += f',y[a+1,{j}]'
        txt_L12_2  += f',y[a,{j}]'
    argu_phase['FCC_A1'] = txt_FCC_A1
    argu_phase['L12_1'] = txt_L12_1
    argu_phase['L12_2'] = txt_L12_2

    for i,input_phase in enumerate(['FCC_A1', 'L12_1', 'L12_2']):
        a = 1
        if i == 2: a = 2

        # cal_f_txtの作成
        if i == 0:
            cal_f_txt += f'    if   i == 0: \n'
            cal_f_txt += f'        f = f_FCC_A1({argu_phase[input_phase]})\n'
        elif i == 1:
            cal_f_txt += f'    elif i == 1  or i == 2 or i == 3 or i == 4: \n'
            cal_f_txt += f'        a = i*2 - 1 \n'
            cal_f_txt += f'        f = f_L12({argu_phase[input_phase]})\n'

        # cal_fy_txtとcal_fyy_txtの作成
        if input_phase == 'FCC_A1':
            cal_fy_txt  += f'    if   a == 0:\n'
            cal_fyy_txt += f'    if   a == 0:\n'
        elif input_phase == 'L12_1':
            cal_fy_txt  += f'    elif a == 1 or a == 3 or a == 5 or a == 7:\n'
            cal_fyy_txt += f'    elif a == 1 or a == 3 or a == 5 or a == 7:\n'
        elif input_phase == 'L12_2':
            cal_fy_txt  += f'    elif a == 2 or a == 4 or a == 6 or a == 8:\n'
            cal_fyy_txt += f'    elif a == 2 or a == 4 or a == 6 or a == 8:\n'
        for j in range(len(solutes)):
            phase_name = input_phase
            if phase_name == 'L12_1': phase_name = 'L12'
            if phase_name == 'L12_2': phase_name = 'L12'
            if j == 0:
                cal_fy_txt  += f'        if   j == {j}: fy = fy_{phase_name}_{solutes[j]}{a}({argu_phase[input_phase]})\n'
                cal_fyy_txt += f'        if   j == {j}: fyy = fyy_{phase_name}_{solutes[j]}{a}({argu_phase[input_phase]})\n'
            else:
                cal_fy_txt  += f'        elif j == {j}: fy = fy_{phase_name}_{solutes[j]}{a}({argu_phase[input_phase]})\n'
                cal_fyy_txt += f'        elif j == {j}: fyy = fyy_{phase_name}_{solutes[j]}{a}({argu_phase[input_phase]})\n'

    cal_f_txt   += f'    return f\n'
    cal_fy_txt  += f'    return fy\n'
    cal_fyy_txt += f'    return fyy\n'

    return cal_f_txt, cal_fy_txt, cal_fyy_txt



# %% [markdown]
# # print_FUNCTION, print_ffyfyy

# %%
def print_FUNCTION(info=False): # Print only the FUCNTIONs to be used｜使用するFUNCTIONのみprintする
    if info:
        print(f'# solvent  : {solvent}')
        print(f'# solutes  : {solutes}')
        print(f'# phases   : {phases}\n')
        if output_style == 'numpy':
            print('import numpy as np')
        elif output_style == 'taichi':
            print('import taichi as ti')

    all_formula = ''
    for phase in phases:
        all_formula += formulate(phase, info=False)
    for i in range(5): # Loop multiple times to find out FUNCTION in FUNCTION｜FUNCTION中FUNCTIONに対応するために複数回ループする
        for FUNCTION_name in list(FUNCTIONS.keys()):
            pattern = rf'(?<=[^\w])({FUNCTION_name})(?=[^\w])' # "not alphabet + FUNCTION_name + not alphabet"｜アルファベット以外+FUNCTION_name+アルファベット以外
            match = re.search(pattern, all_formula)
            if match:
                all_formula += FUNCTIONS[FUNCTION_name]
    for FUNCTION_name in list(FUNCTIONS.keys()):
        pattern = rf'(?<=[^\w])({FUNCTION_name})(?=[^\w])' # アルファベット以外+FUNCTION_name+アルファベット以外
        match = re.search(pattern, all_formula)
        if match:
            function_text = FUNCTIONS[FUNCTION_name]
            if output_style == 'numpy':
                function_text = function_text.replace('ti.log','np.log').replace('@ti.func','')
            print(function_text)  
    
def print_ffyfyy(info=False):
    for phase in phases:
        phase_formula = formulate(phase, info=info)
        f_txt, fy_txt, fyy_txt = construct_f_fy_fyy(phase_formula, phase, info=info)
        if output_style == 'numpy':
            f_txt = f_txt.replace('ti.log','np.log').replace('@ti.func','').replace(':int','').replace(' -> ti.f64','').replace(':ti.f64','')
            fy_txt = fy_txt.replace('ti.log','np.log').replace('@ti.func','').replace(':int','').replace(' -> ti.f64','').replace(':ti.f64','')
            fyy_txt = fyy_txt.replace('ti.log','np.log').replace('@ti.func','').replace(':int','').replace(' -> ti.f64','').replace(':ti.f64','')
        print(f_txt)
        print(fy_txt)
        print(fyy_txt)
        
    if useful_function_for_PFM:
        cal_f_txt, cal_fy_txt, cal_fyy_txt = construct_function_for_PFM()
        if userful_function_for_PFM_Ni_superalloy:
            cal_f_txt, cal_fy_txt, cal_fyy_txt = construct_function_for_PFM_superalloy()
        if output_style == 'numpy':
            cal_f_txt   = cal_f_txt.replace('@ti.func','').replace(':int','').replace(' -> ti.f64','').replace(':ti.f64','')
            cal_fy_txt  = cal_fy_txt.replace('@ti.func','').replace(':int','').replace(' -> ti.f64','').replace(':ti.f64','')
            cal_fyy_txt = cal_fyy_txt.replace('@ti.func','').replace(':int','').replace(' -> ti.f64','').replace(':ti.f64','')
        print(cal_f_txt)
        print(cal_fy_txt)
        print(cal_fyy_txt)
    
    if full_TDB_text:
        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                print(f'# {line.strip()}')
                


# %% [markdown]
# # print text

# %%
print_FUNCTION(info=True)
print_ffyfyy(info=supplementary_infomation)
