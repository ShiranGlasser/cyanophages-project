from record.record import fields


def get_main_attributes(records):
    keys = sorted(fields)
    main_attributes = {}
    for record in records:
        main_attributes[record.record_id] = get_attributes(keys, record)
    return main_attributes


def get_attributes(keys, record):
    attributes = {}
    attributes_keys = record.main_attributes.keys()
    for key in keys:
        if key in attributes_keys:
            attributes[key] = record.main_attributes[key]
        else:
            attributes[key] = ""
    return attributes

def stats_amino_acids(codon_table):
    dic = {}
    for elemnt in codon_table:
        if codon_table[elemnt] in dic.keys():
            dic[codon_table[elemnt]] = dic[codon_table[elemnt]] + 1
        elif codon_table[elemnt] == '_':
            if 'END' in dic.keys():
                dic['END'] = dic['END'] + 1
            else:
                dic['END'] = 1
        else:
            dic[codon_table[elemnt]] = 1
    return dic