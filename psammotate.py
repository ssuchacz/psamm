import csv
from collections import defaultdict, OrderedDict
from psamm.datasource.native import NativeModel, ReactionEntry
import argparse
import logging
import re
from psamm.expression import boolean
from psamm_import import importer


def run():
    logging.basicConfig()
    parser = argparse.ArgumentParser(description='File Path to the model.yaml file')
    parser.add_argument('--model', help='path to model.yaml file', metavar="FILE")
    parser.add_argument('--pan_res', help='path to appearance_2 file from panscript', metavar="FILE")
    args = parser.parse_args()
    mm_file = args.model
    pan_result = args.pan_res
    print(mm_file)
    tdict = app_reader(pan_result, 0, 1)
    print(tdict)
    trans_genes_dict = model_loader(mm_file, tdict)
    model_export(mm_file, trans_genes_dict)


def app_reader(app_file, query, template):
    '''This function will read in the app file and produce the mapping dictionary

    The input is the app_file argument that the user enters, the query genome, which
    is the number of the genome that the user wishes to make a model of, and the
    number of the template genome which the user wishes to use to create the new model
    from.
    '''

    app_file_o = open(app_file, mode='rU')
    new_list = {}
    trans_dict = defaultdict(list)
    transl_dict = {}
    for x, row in enumerate(csv.reader(app_file_o, delimiter='\t')):
        temp_l = []
        quer_l = []
        if x == 0:
            new_list['template'] = (row[template])
            new_list['query'] = (row[query])
        else:
            temp = row[template]
            quer = row[query]
            temp_split = []
            if ',' in temp:
                temp_split = temp.split(',')
            else:
                temp_split.append(temp)
            for i in temp_split:
                if ':' in i:
                    i = i.split(':')[1]
                else:
                    i = i
                temp_l.append(i)
            if ',' in quer:
                quer_split = quer.split(',')
            else:
                quer_split = [quer]
            for i in quer_split:
                if ':' in i:
                    i = i.split(':')[1]
                else:
                    i = i
                quer_l.append(i)
        for i in temp_l:
            i = i.strip()
            for j in quer_l:
                j = j.strip()
                trans_dict[i].append(j)
        for key, value in trans_dict.iteritems():
            value = str(value)
            value = value.strip('[')
            value = value.strip(']')
            value = value.strip()
            value = value.strip(',')
            value = value.replace(',', ' or')
            transl_dict[key] = value
    return trans_dict

def model_loader(model_file, translation_dict):
    new_translation_dict = {}
    translated_genes = {}
    model = NativeModel.load_model_from_path(model_file)
    print(model.name)
    target_genes_l = {}
    target_model_reactions = []
    translation_dict.pop('-', None)
    for key, value in translation_dict.iteritems():
        value_s = None
        for i in value:
            if value_s is None:
                value_s = i
            else:
                value_s = value_s + ' or ' + i
        if len(value) >=2:
            value_s = '(' + value_s + ')'
        new_translation_dict[key] = value_s
        new_translation_dict['Sputw3181_3559'] = '-'
        new_translation_dict['TM1418'] = '-'
        new_translation_dict['TM1438'] = '-'
    for value in translation_dict.values():
        for i in value:
            target_genes_l[i] = True
    target_genes_l.pop('-', None)
    target_genes_l.pop('\'-\'', None)
    target_genes_l['\'-\''] = False
    target_genes_l['-'] = False
    target_genes_l['S_0001'] = True
    for entry in model.parse_reactions():
        print(entry.genes)
        if entry.genes is None:
            target_model_reactions.append(entry.id)
            translated_genes[entry.id] = None
        else:
            genes = entry.genes
            genes = genes.strip('_i')
            ###Uncoment line relative to target genome
            ###Needed for modification of the gene IDs so that the ong nomenclature matches the kegg nomenclature
            #genes = genes.replace('SO','SO_')
            #genes = genes.replace('MR4', 'Shewmr4_')
            #genes = genes.replace('W3181_', 'Sputw3181_')
            genes = genes.replace('TM_', 'TM')
            genes_1 = entry.genes
            for key, value in new_translation_dict.iteritems():
                #print(genes, key, value)
                genes = re.sub(key, value, genes)
            print('###')
            print(genes_1)
            print(genes)
            e = boolean.Expression(genes)

            e_1 = e.substitute(lambda v: target_genes_l.get(v.symbol, v))
            #if e_1.value is True:
            #    target_model_reactions.append(entry.id)
            #    translated_genes[entry.id] = genes

            #or
            #translated_genes[entry.id] = [genes_1, genes, e_1.value]
            translated_genes[entry.id] = [genes_1, genes, False]


    return translated_genes


def model_export(model_file, trans_genes):
    f = open('Target_genome_mapping.tsv', mode='w')
    for key, value in trans_genes.iteritems():
        if value is not None:
            f.write('{}\t{}\t{}\t{}\n'.format(key, value[0], value[1], value[2]))
        else:
            f.write('{}\t{}\t{}\t{}\n'.format(key, 'None', 'None', 'None'))




run()
