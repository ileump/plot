#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Imports
import os, sys, datetime
import logging
import argparse, optparse
import numpy as np
import pandas as pd
import pyparsing as pp
import pathlib
import typing

## lynxx change cwd
# cwd = os.getcwd()
# from lynx.controllers.converter import Converter
# from lynx.utils.log import app_logger
# os.chdir(cwd)

def fa_parser() -> pp.ParserElement:
    """
    Parse fatty acyl with finer resolution.
    """
    # fatty acyl rules
    carbon = pp.Word(pp.nums)('carbon')
    db = pp.Word(pp.nums)('db_count')

    hydroxyl = pp.oneOf('h d t')('hydroxyl')
    ether = pp.oneOf('p e')('ether')

    lcb_core = pp.Combine(hydroxyl + carbon + ':' + db)
    lcb = lcb_core('lcb')

    fa_core = pp.Combine(carbon +':' + db + pp.Optional(ether) + pp.Optional(hydroxyl))
    fa = fa_core | lcb

    # fa.runTests("""
    # 32:2
    # d18:1
    # 22:0h
    # 40:6p
    # 32:2e
    # """)
    return(fa)

def sort_fa(x: str):
    ## assume carbon has 2 digits at most and db has 2 digit at most
    x_split = x.split(':')
    ## handle prefix O-, P-
    x0_split = x_split[0].split('-')

    ## with prefix
    if len(x0_split) > 1:
        res = x0_split[0] + '-' + x0_split[1].rjust(2, '0') + x_split[1].rjust(2, '0')
    else:
        ## use Z- as prefix so that O- and P- comes first
        res = 'Z-' + x_split[0].rjust(2, '0') + x_split[1].rjust(2, '0')

    return(res)

def fa_to_M0(fa: str) -> str:
    faparser = fa_parser()
    tokens = faparser.parseString(fa)
    if 'lcb' in tokens.keys():
        if tokens['lcb']['hydroxyl'] == 'd':
            res = '{carbon}:{db_count};O2'.format(
                carbon = tokens['lcb']['carbon'],
                db_count = tokens['lcb']['db_count']
            )
        elif tokens['lcb']['hydroxyl'] == 't':
            res = '{carbon}:{db_count};O3'.format(
                carbon = tokens['lcb']['carbon'],
                db_count = tokens['lcb']['db_count']
            )
    elif 'hydroxyl' in tokens.keys():
        if tokens['hydroxyl'] == 'h':
            res = '{carbon}:{db_count};O'.format(
                    carbon = tokens['carbon'],
                    db_count = tokens['db_count']
                )
    else:
        res = '{carbon}:{db_count}'.format(
            carbon = tokens['carbon'],
            db_count = tokens['db_count']
        )

    if fa[-1] == 'p':
        res = 'P-' + res
    elif fa[-1] == 'e':
        res = 'O-' + res

    return(res)

def fas_to_B0(*fas: str) -> str:
    faparser = fa_parser()

    n_carbon = 0
    n_db = 0
    n_O_suffix = 0
    prefix = ''

    for fa in fas:
        tokens = faparser.parseString(fa)
        if 'lcb' in tokens.keys():
            if tokens['lcb']['hydroxyl'] == 'd':
                n_carbon += int(tokens['lcb']['carbon'])
                n_db += int(tokens['lcb']['db_count'])
                n_O_suffix += 2
            elif tokens['lcb']['hydroxyl'] == 't':
                n_carbon += int(tokens['lcb']['carbon'])
                n_db += int(tokens['lcb']['db_count'])
                n_O_suffix += 3
        elif 'hydroxyl' in tokens.keys():
            if tokens['hydroxyl'] == 'h':
                n_carbon += int(tokens['carbon'])
                n_db += int(tokens['db_count'])
                n_O_suffix += 1
        else:
            n_carbon += int(tokens['carbon'])
            n_db += int(tokens['db_count'])

        if fa[-1] == 'p':
            prefix = 'O-'
            n_db += 1
        elif fa[-1] == 'e':
            prefix = 'O-'

    if n_O_suffix > 1:
        res = '{carbon}:{db_count};O{n_O}'.format(
            carbon = n_carbon,
            db_count = n_db,
            n_O = n_O_suffix
        )
    elif n_O_suffix == 1:
        res = '{carbon}:{db_count};O'.format(
            carbon = n_carbon,
            db_count = n_db
        )
    else:
        res = '{carbon}:{db_count}'.format(
            carbon = n_carbon,
            db_count = n_db
        )
    
    res = prefix + res
    return(res)

def init_lipidall_parser() -> pp.ParserElement:
    
    # fatty acyl rules
    carbon = pp.Word(pp.nums)
    db = pp.Word(pp.nums)

    hydroxyl = pp.oneOf('h d t')
    ether = pp.oneOf('p e')

    lcb_core = pp.Combine(hydroxyl + carbon + ':' + db)
    lcb = lcb_core

    fa_core = pp.Combine(carbon +':' + db + pp.Optional(ether) + pp.Optional(hydroxyl))
    fa = fa_core | lcb
    fa_sum = fa_core

    fa2_unsorted = fa('fa1') + pp.Suppress('/') + fa('fa2')
    fa2 = fa2_unsorted

    # FFA
    ffa = pp.Literal('FFA')('hg') + fa('fa1')

    # ffa.runTests("""
    # FFA22:6
    # """)

    # carnitine
    car = fa('fa1') + '-' + pp.Literal('carnitine')('hg')

    # car.runTests("""
    # 18:0-carnitine
    # """)

    # glycerolipid rules
    gl_two = pp.Literal('DAG')('hg') + fa_sum('fa_sum') + pp.Suppress('(') + fa2 + pp.Suppress(')') 
    gl_three = pp.Literal('TAG')('hg') + fa_sum('fa_sum') + pp.Suppress('(') + fa('fa1_3') + pp.Suppress(')')
    gl = gl_two | gl_three

    # gl.runTests("""
    # TAG42:0(14:0)
    # DAG40:5(18:0/22:5)
    # """)

    # phospholipid rules 
    pl_one_hg = pp.Literal('LPA') | 'LPC' | 'LPE' | 'LPG' | 'LPI' | 'LPS' | 'LysoPC'
    pl_one = pl_one_hg('hg') + fa('fa1')

    pl_two_hg = pp.Literal('PA') | 'PC' | 'PE' | 'PG' | 'PI' | 'PS' | 'BMP'
    pl_two = pl_two_hg('hg') + fa_sum('fa_sum') + pp.Optional(pp.Suppress('(') + fa2 + pp.Suppress(')'))

    pl_four_hg = pp.Literal('CL')
    pl_four = pl_four_hg('hg') + fa_sum('fa_sum') + pp.Suppress('(') + fa('fa1_4') + pp.Suppress(')')

    pl = pl_one | pl_two | pl_four

    # pl.runTests("""
    # PG32:2
    # BMP32:2
    # CL66:4(16:1)
    # PE32:0
    # PI 34:2
    # PC32:1(16:0/16:1)
    # LPI16:0
    # LysoPC16:0
    # LPE18:1p
    # PC38:4p(18:0/20:4)
    # PC32:2e
    # """)

    # sphingolipid rules
    ganglioside = pp.Literal('Gb3') | 'GM3' | 'GM2' | 'GM1'
    sl_hg  = pp.Literal('Cer') | 'SM' | 'GluCer' | 'LacCer' | 'SL' | ganglioside

    sl_lcb_subspecies = lcb('lcb') + '/' + fa('fa1')
    sl_lcb = sl_lcb_subspecies

    sl_base = pp.Literal('Sph')

    sl = sl_hg('hg') + sl_lcb | sl_base('hg') + lcb('lcb')

    # sl.runTests("""
    # GM3 d18:1/16:0
    # SM d18:1/14:0
    # Cer d18:1/16:0
    # GluCer d18:1/16:0
    # LacCer d18:1/14:0
    # Gb3 d18:1/16:0
    # SL d18:1/22:0h
    # Sph d18:1
    # """)

    # sterol rules
    st_sub1_hg = pp.Literal('CE')
    st_sub1 = st_sub1_hg('hg') + pp.Suppress('-') + fa('fa1')

    st = st_sub1

    # st.runTests("""
    # CE-22:0
    # """)

    ## all classes
    lipid = ffa | car | gl | pl | sl | st

    return(lipid)

def parseAction_factory() -> typing.Callable:
    """
    """
    def lipidall_lipidlynxx_parseAction(tokens: dict) -> dict:
        class_mapping = {
            'LysoPA': 'LPA',
            'LysoPC': 'LPC',
            'LysoPE': 'LPE',
            'LysoPG': 'LPG',
            'LysoPI': 'LPI',
            'LysoPS': 'LPS',
            'BMP': 'LBPA',
            'SL': 'SHexCer',
            'carnitine': 'CAR',
            'Sph': 'SPB',
            'GluCer': 'GlcCer',
            'DAG': 'DG',
            'TAG': 'TG',
            'FFA': 'FA'
        }

        lipid_class = tokens['hg']
        
        if lipid_class in class_mapping.keys():
            lipid_class = class_mapping[tokens['hg']]

        if 'lcb' in tokens.keys():
            if 'fa1' in tokens.keys():
                tokens['lipidlynxx_M0'] =  '{lipid_class}({lcb}_{fa1})'.format(
                    lipid_class = lipid_class,
                    lcb = fa_to_M0(tokens['lcb']),
                    fa1 = fa_to_M0(tokens['fa1'])
                )
                tokens['lipidlynxx_B0'] = '{lipid_class}({fa_B0})'.format(
                    lipid_class = lipid_class,
                    fa_B0 = fas_to_B0(tokens['lcb'], tokens['fa1'])
                )
            else:
                tokens['lipidlynxx_M0'] =  '{lipid_class}({lcb})'.format(
                    lipid_class = lipid_class,
                    lcb = fa_to_M0(tokens['lcb'])
                )
                tokens['lipidlynxx_B0'] = '{lipid_class}({fa_B0})'.format(
                    lipid_class = lipid_class,
                    fa_B0 = fas_to_B0(tokens['lcb'])
                )
        elif 'fa_sum' in tokens.keys() and \
            len(set(['fa1', 'fa2', 'fa3', 'fa4']) & set(tokens.keys())) > 0 and \
            tokens['fa_sum'][-1] in ['e', 'p']:
                tokens['lipidlynxx_B0'] = '{lipid_class}({fa})'.format(
                    lipid_class = lipid_class,
                    fa = fas_to_B0(tokens['fa_sum'])
                )
                tokens['fa1_2'] = tokens['fa1']
                tokens['fa2_2'] = tokens['fa2']
        elif len(set(['fa1', 'fa2', 'fa3', 'fa4']) & set(tokens.keys())) > 0:
            fas = []
            fas_M0 = []
            for k in ['fa1', 'fa2', 'fa3', 'fa4']:
                if k in tokens.keys():
                    fas.append(tokens[k])
                    fas_M0.append(fa_to_M0(tokens[k]))

            fas_M0 = list(filter(lambda x: x != '0:0', fas_M0))
            fas_M0.sort(key = sort_fa)
            tokens['lipidlynxx_M0'] =  '{lipid_class}({fa})'.format(
                lipid_class = lipid_class,
                fa = '_'.join(fas_M0)
            )
            tokens['lipidlynxx_B0'] = '{lipid_class}({fa})'.format(
                lipid_class = lipid_class,
                fa = fas_to_B0(*fas)
            )
        elif 'fa_sum' in tokens.keys():
            tokens['lipidlynxx_B0'] = '{lipid_class}({fa})'.format(
                lipid_class = lipid_class,
                fa = fas_to_B0(tokens['fa_sum'])
            )

    return(lipidall_lipidlynxx_parseAction)

def convert_lipidname(lipidall_names: pd.Series) -> pd.DataFrame:
    """
    """
    lipidall_parser = init_lipidall_parser()
    parser = lipidall_parser.addParseAction(parseAction_factory())

    # parser.runTests("""
    # PG32:2
    # BMP32:2
    # SL d18:1/22:0
    # SL d18:1/22:0h
    # CL66:4(16:1)
    # FFA22:6
    # PE40:6p
    # PE32:1
    # PI 34:2
    # GM3 d18:1/16:0
    # PC32:2(16:1/16:1)
    # PC36:4p(16:0/20:4)
    # LPI16:0
    # LPE18:1p
    # LPA18:2
    # PA32:2
    # LPS16:0
    # PS 34:2
    # LysoPC16:0
    # PC32:2e
    # PC34:2p
    # PC40:5
    # SM d18:1/14:1
    # 14:0-carnitine
    # Sph d18:1
    # Cer d18:1/16:0
    # GluCer d18:1/16:0
    # LacCer d18:1/14:0
    # Gb3 d18:1/16:0
    # DAG28:0(14:0/14:0)
    # TAG42:0(14:0)
    # CE-16:1
    # """)

    inputfile = "F:/Lipidall/projects/Reports/2020/2020-155-B-01 人类肝癌膀胱癌/2020-155-B-01 I 人类肝癌脂质组/data/20210311/2020-155-B-01-I.csv"
    inputfile = "F:/Lipidall/projects/Reports/2020/2020-155-B-01 人类肝癌膀胱癌/2020-155-B-01-III 人类膀胱癌脂质组/data/20200914/2020-155-B-01-III.csv"
    data = pd.read_csv(inputfile)

    ## rename first column
    data.rename(columns = {data.columns[0]:"rownames"}, inplace=True)
    
    data_lipidname = pd.DataFrame(
        [('', '', '', '', '') for _ in range(data.shape[0])],
        columns = ['lipidall', 'lipidlynxx_B0', 'lipidlynxx_M0', 'lipidmaps_abbrev', 'biopan']
    )

    for i in range(len(data['rownames'])):
        data_lipidname.iloc[i]['lipidall'] = data.iloc[i]['rownames']
        try:
            out = parser.parseString(data.iloc[i]['rownames'])
            if 'lipidlynxx_B0' in out.keys():
                data_lipidname.iloc[i]['lipidlynxx_B0'] =  out['lipidlynxx_B0']
            if 'lipidlynxx_M0' in out.keys():
                data_lipidname.iloc[i]['lipidlynxx_M0'] =  out['lipidlynxx_M0']
            if 'lipidmaps_abbrev' in out.keys():
                data_lipidname.iloc[i]['lipidmaps_abbrev'] =  out['lipidmaps_abbrev']
            if 'biopan' in out.keys():
                data_lipidname.iloc[i]['biopan'] =  out['biopan']
        except pp.ParseException:
            out = ''

    data_lipidname.to_csv('output/2020-155-B-01-I-mapping.csv')

# Module Functions and Classes
def main(argv=None):
    """Main script function.
    """

    # Parse command line arguments
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", help="input file")
    parser.add_argument("outputfile", help="input file")
    args = parser.parse_args()

    ## current directory is modified in some function
    inputfile = pathlib.Path.absolute(pathlib.Path(args.inputfile))
    outputfile = pathlib.Path.absolute(pathlib.Path(args.outputfile))

    print(outputfile)

    if (inputfile.exists()):
        print("Message: input file ", inputfile, 'exists.')
    else:
        print("Message: input file ", inputfile, 'does not exists.')
        raise Exception('Input file is not found.')

    # logging
    logging.basicConfig(
        filename=pathlib.Path.joinpath(
            pathlib.Path("."),
            pathlib.Path("lipid_convert.log")),
        level = logging.INFO
    )
    logging.info("{:%Y-%m-%d %H:%M:%S}".format(datetime.datetime.now()) + ' ' +
                 'Started\n')

    data = pd.read_csv(inputfile)

    ## rename first column
    data.rename(columns = {data.columns[0]:"rownames"}, inplace=True)
    
    data_lipidname = pd.DataFrame(
        [('', '', '', '', '', '', '') for _ in range(data.shape[0])],
        columns = ['lipidall', 'lipidlynxx_B0', 'lipidlynxx_M0', 'fa1_2', 'fa2_2', 'fa1_3', 'fa1_4']
    )
    
    lipidall_parser = init_lipidall_parser()
    parser = lipidall_parser.addParseAction(parseAction_factory())

    for i in range(data.shape[0]):
        data_lipidname.iloc[i]['lipidall'] = data.iloc[i]['rownames']
        try:
            tokens = parser.parseString(data.iloc[i]['rownames'])
            if 'lipidlynxx_B0' in tokens.keys():
                data_lipidname.iloc[i]['lipidlynxx_B0'] =  tokens['lipidlynxx_B0']
            if 'lipidlynxx_M0' in tokens.keys():
                data_lipidname.iloc[i]['lipidlynxx_M0'] =  tokens['lipidlynxx_M0']
            if 'fa1_2' in tokens.keys():
                data_lipidname.iloc[i]['fa1_2'] = tokens['fa1_2']
            if 'fa2_2' in tokens.keys():
                data_lipidname.iloc[i]['fa2_2'] = tokens['fa2_2']  
            if 'fa1_3' in tokens.keys():
                data_lipidname.iloc[i]['fa1_3'] = tokens['fa1_3']
            if 'fa1_4' in tokens.keys():
                data_lipidname.iloc[i]['fa1_4'] = tokens['fa1_4']
        except pp.ParseException:
            pass

    data_lipidname.to_csv(outputfile)

# Check to see if this file is the "__main__" script being executed
if __name__ == '__main__':
    main()
