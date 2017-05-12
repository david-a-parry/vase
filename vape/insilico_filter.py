

class InsilicoFilter(object):
    ''' Stores in silico prediction formats for VEP and assesses with a
        score for a particular program indicates pathogenicity or not.
        Data on in silico formats recognised are stored in
        ../data/vep_insilico_pred.tsv
    '''
    
    
    def __init__(self, programs=None):
        ''' Initialize with a list of program names to use as filters.
            Optionally specify score criteria for filtering as in the
            the following examples:
                                    FATHMM_pred=D
                                    MutationTaster_pred=A
                                    MetaSVM_rankscore=0.8
        '''
        self.pred_filters = {}
        self.score_filters = {}
        default_progs = {}
        case_insensitive = {}
        lower_more_damaging = set()
        with open ("../data/vep_insilico_pred.tsv", encoding='UTF-8') as insilico_data:
            for line in insilico_data:
                if line.startswith('#'):
                    continue
                cols = line.rstrip().split('\t')
                case_insensitive[cols[0].lower()] = cols[0]
                if cols[0] in default_progs:
                    if cols[2] in default_progs[cols[0]]:
                        default_progs[cols[0]][cols[2]].append(cols[1])
                    elif cols[2] != 'score' :
                        default_progs[cols[0]][cols[2]] = [cols[1]]
                    else:
                        raise Exception("Error in data/vep_insilico_pred.tsv: should" +
                                        " only have one entry for score prediction " + 
                                        "'{}'" .format(cols[0]))
                else:
                    if cols[2] != 'score':
                        default_progs[cols[0]] = {'type' : 'pred'}
                        default_progs[cols[0]][cols[2]] = [cols[1]]
                    else:
                        default_progs[cols[0]] = {'type' : 'score', 
                                                  'default' : float(cols[1])}
                if len(cols) >= 4:
                    if cols[3] == 'lower=damaging':
                        lower_more_damaging.add(cols[0])

        for prog in programs:
            split = prog.split('=')
            pred = None
            if len(split) > 1:
                prog = split[0]
                pred = split[1]
            if prog.lower() in case_insensitive:
                prog = case_insensitive[prog.lower()]
            else:
                raise Exception("ERROR: in silico prediction program '{}' "
                                .format(prog) + "not recognised.")
            if pred is not None:
                if default_progs[prog]['type'] == 'score':
                    try:
                        score = float(pred)
                        self.score_filters[prog] = score
                    except ValueError:
                        raise Exception("ERROR: {} score must be numeric. " 
                                        .format(prog) + "Could not convert " +
                                        "value '{}' to a number.".format(pred))
                elif (pred in default_progs[prog]['default'] or 
                     pred in default_progs[prog]['valid']):
                    if prog in self.pred_filters[prog]:
                        self.pred_filters[prog].add(pred)
                    else:
                        self.pred_filters[prog] = set(pred)
                else:
                    raise Exception("ERROR: score '{}' not " .format(pred) +
                                    "recognised as valid for in silico " + 
                                    "prediction program '{}' ".format(prog))
                        

    def check_pred(self, prog, pred):
        ''' 
            Return True if prediction matches filter for given program,
            otherwise return False.
        
            Args:
                prog: Name of in silico program - must have been 
                      specified when initialized
                
                pred: Prediction value to check.

        '''
        pass
