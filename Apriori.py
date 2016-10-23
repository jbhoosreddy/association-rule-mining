from __future__ import division
import itertools
import pickle
import re

# from multiprocessing import Pool
# from pathos.multiprocessing import ProcessingPool as Pool


template1 = ['RULE HAS ANY OF G6_UP', 'RULE HAS 1 OF G1_UP', 'RULE HAS 1 OF (G1_UP, G10_DOWN)', 'BODY HAS ANY OF G6_UP',
             'BODY HAS NONE OF G72_UP', 'BODY HAS 1 OF (G1_UP, G10_DOWN)', 'HEAD HAS ANY OF G6_UP',
             'HEAD HAS NONE OF (G1_UP, G6_UP)', 'HEAD HAS 1 OF (G6_UP, G8_UP)', 'RULE HAS 1 OF (G1_UP, G6_UP, G72_UP)',
             'RULE HAS ANY OF (G1_UP, G6_UP, G72_UP)']

template2 = ['SIZE OF RULE >= 3', 'SIZE OF BODY >= 2', 'SIZE OF HEAD >= 2']

template3 = ['BODY HAS ANY OF G1_UP AND HEAD HAS 1 OF G59_UP', 'BODY HAS ANY OF G1_UP OR HEAD HAS 1 OF G6_UP',
             'BODY HAS 1 OF G1_UP OR HEAD HAS 2 OF G6_UP', 'HEAD HAS 1 OF G1_UP AND BODY HAS 0 OF DISEASE',
             'HEAD HAS 1 OF DISEASE OR RULE HAS 1 OF (G72_UP, G96_DOWN)',
             'BODY HAS 1 of (G59_UP, G96_DOWN) AND SIZE OF RULE >= 3']


def load_data(path):

    file_handle = open(path)
    data = dict()
    unique = list();
    for line in file_handle:
        line = line.rstrip()
        line = line.split('\t')
        symptoms = list()
        x = 1
        for gene in line[1:-1]:
            gene = 'G' + str(x) + '_' + gene.upper()
            symptoms.append(gene)
            x += 1

        symptoms.append(line[-1])

        data[line[0]] = symptoms
        unique.extend(symptoms)
    unique = set(unique)
    file_handle.close()
    return data, unique, len(data.keys())


def create_structure(data, genes):
    struct = dict()
    for gene in genes:
        for key, value in data.items():
            if gene in value:
                if gene not in struct.keys():
                    struct[gene] = list()
                struct[gene].append(key)
    return struct


def prune(itemsets, total, support, length):
    itemset = dict(itemsets[length]['data'])
    pruned = list()
    if len(itemset) == 0:
        return itemsets
    minThreshold = support * total
    for k, v in itemset.items():
        if len(v) < minThreshold:
            del itemset[k]
            pruned.append(k)
    itemsets[length]['data'] = itemset
    itemsets[length]['pruned'] = pruned
    return itemsets


def populate(itemsets, length):
    original = itemsets[0]['data']
    prev = itemsets[length-1]['data']
    if not len(prev):
        return itemsets
    prevPruned = itemsets[length-1]['pruned']
    temp = dict()
    combinations = create_combinations(prev.keys(), prevPruned, length+1)
    for subset in combinations:
        count = 1
        intersect_set = None
        for key in subset:
            if count == 1:
                intersect_set = set(original[key])
                count += 1
            else:
                intersect_set = intersect_set.intersection(set(original[key]))
                if not len(intersect_set):
                    break
        if len(intersect_set):
            temp[','.join(subset)] = intersect_set

    itemset = dict()
    itemset['data'] = temp
    itemsets.append(itemset)
    return itemsets


def get_sample(itemset, key):
    if itemset is None:
        return
    if key in itemset.keys():
        return itemset[key]
    return None


def create_combinations(aKeys, aPruned, length):
    log = False
    if length == 3:
        log = True

    keys = map(lambda x: x.split(',') if ',' in x else [x], aKeys)

    pruned = map(lambda x: x.split(',') if ',' in x else [x], aPruned)

    keys = [item for sublist in keys for item in sublist]
    keys = set(keys)
    combinations = itertools.combinations(keys, length)
    keys = filter(lambda x: not is_pruned(x, pruned), combinations)
    return keys


def is_pruned(key, pruned):
    flag = False
    for aPrune in pruned:
        flag = all(p in key for p in aPrune)
        if flag:
            return flag
    return flag




def save(itemsets):
    with open('itemsets.pickle', 'wb') as handle:
        pickle.dump(itemsets, handle)


def get(itemsets, length=None, key=None):
    if length is None:
        return itemsets
    if key is None:
        return itemsets[length]['data']
    return itemsets[length]['data'][key]


def generate_rule(body, keys, id):
    body = set(body)
    keys = set(keys)
    head = keys - keys.intersection(body)
    return create_rule(list(body), list(head), id)


def create_rule(body, head, id):
    rule = dict()
    rule['BODY'] = body
    rule['HEAD'] = head
    rule['RULE'] = list()
    rule['RULE'].extend(head)
    rule['RULE'].extend(body)
    rule['RULE'] = set(rule['RULE'])
    rule['ID'] = id
    return rule


def classify_rules_from_query_debug(query, all_rules):
    return classify_rules_from_query(query, all_rules)


def classify_rules_from_query(query, all_rules):
    if 'AND' not in query and 'OR' not in query:
        return generate_rules_from_query(query, all_rules)
    if 'AND' in query:
        method = 'intersection'
        query = map(lambda x: x.strip(), query.split('AND'))
    elif 'OR' in query:
        method = 'union'
        query = map(lambda x: x.strip(), query.split('OR'))
    rules = map(lambda q: generate_rules_from_query(q, all_rules), query)
    return join(rules[0], rules[1], method, all_rules)


def generate_rules_from_query(query, all_rules):
    if 'HAS' in query:
        return rules_from_template1(query, all_rules)
    elif 'SIZE OF' in query:
        return rules_from_template2(query, all_rules)


def join(rules1, rules2, method, all_rules):
    idx1 = set(map(lambda x: x['ID'], rules1))
    # print rules1, idx1
    idx2 = set(map(lambda x: x['ID'], rules2))
    if method == 'union':
        idx = idx1.union(idx2)
    elif method == 'intersection':
        idx = idx1.intersection(idx2)
    return filter(lambda x: x['ID'] in idx, all_rules)


def validate_rule(rule, items, portion, constraint):
    if constraint == 'NONE':
        if not len(filter(lambda item: item in rule[portion], items)):
            return True
        return False
    if constraint == 'ANY':
        if len(filter(lambda item: item in rule[portion], items)) > 0:
            return True
        return False
    else:
        if len(filter(lambda item: item in rule[portion], items)) == int(constraint):
            return True
        return False


def rules_from_template1(query, all_rules):
    items = re.search("\(.*\)", query)
    query = query.split(' ')
    first = query[0]
    number = query[2]
    if items:
        items = map(lambda x: x.strip(), items.group()[1:-1].split(','))
    else:
        items = [query[-1]]

    return filter(lambda rule: validate_rule(rule, items, first, number), all_rules)


def rules_from_template2(query, all_rules):
    query = query.split(' ')
    number = int(query.pop())
    portion = query[-2]
    return filter(lambda rule: len(rule[portion]) >= number, all_rules)


def generate_rules(itemsets, total, confidence):
    rules = list()
    # confidence = confidence * total
    # print 'In Rule generation'
    max_length = len(itemsets)
    id = 0
    for size in range(max_length-1, 1, -1):
        # print size
        itemset = get(itemsets, size-1)
        keys = itemset.keys()
        # print keys
        for aKey in keys:
            i = 0
            total_numerator = len(itemset[aKey])
            # print 'total_numerator', total_numerator
            while size-i > 1:
                i += 1
                key = aKey.split(',')
                combination = itertools.combinations(key, size-i)
                den_itemset = get(itemsets, size-i-1)
                for subset in combination:
                    if len(subset) > 1:
                        sub_permutation = itertools.permutations(subset, size-i)
                        for sub_subset in sub_permutation:
                            subset_key = ','.join(sub_subset)
                            if subset_key in den_itemset.keys():
                                den = den_itemset[subset_key]
                                break
                    else:
                        den = den_itemset[subset[0]]

                    total_den = len(den)
                    # print total_den
                    rule_confidence = total_numerator/total_den
                    if rule_confidence >= confidence:
                        rule = generate_rule(subset, key, id)
                        id += 1
                        rule['CONFIDENCE'] = rule_confidence
                        rules.append(rule)
    return rules


def print_rules(rules):
    for rule in rules:
        print rule['BODY'], '-->', rule['HEAD'], '(', rule['CONFIDENCE'], ')', '[', rule['ID'], ']'


def main(support, confidence, cached=False):
    data, genes, total = load_data('/Users/vinaygoyal/Desktop/gene_expression.txt')
    if cached:
        with open('itemsets.pickle', 'rb') as handle:
            itemsets = pickle.load(handle)
    else:
        itemsets = list()
        itemset = dict()
        itemset['data'] = create_structure(data, genes)
        itemsets.append(itemset)
        for i in range(5):
            # print 'Original', i, len(itemsets[i]['data'])
            itemsets = prune(itemsets, total, support, i)
            print 'Prune', i, len(itemsets[i]['data'])
            if i == 2:
                print 'length 3 data', (itemsets[i]['data'])
            itemsets = populate(itemsets, i+1)
            print 'Populate', (i+1), len(itemsets[i+1]['data'])

            if not len(itemsets[i+1]['data']):
                save(itemsets)
                break
    all_rules = generate_rules(itemsets, total, confidence)
    # print_rules(all_rules)
    for query in template1:
        print "Query:", query
        rules = classify_rules_from_query(query, all_rules)
        print 'Rule length:', len(rules)
        # print_rules(rules)
    for query in template2:
        print "Query:", query
        rules = classify_rules_from_query(query, all_rules)
        # print_rules(rules)
        print 'Rule length:', len(rules)
    for query in template3:
        print "Query:", query
        rules = classify_rules_from_query(query, all_rules)
        print 'Rule length:', len(rules)
        # print_rules(rules)
    debug = False
    while True:
        print
        query = raw_input("Enter query: ")
        if query == "DEBUG":
            debug = not debug
            continue
        if not query: break
        if debug:
            rules = classify_rules_from_query_debug(query, all_rules)
        else:
            rules = classify_rules_from_query(query, all_rules)
        print 'Rule length:', len(rules)
        # print_rules(rules)

# Main function --> main(support , confidence, cache)
main(0.3, 0.6, False)
