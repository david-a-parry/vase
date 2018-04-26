import requests
import time
import logging

server = "https://rest.ensembl.org"
grch37_server = "http://grch37.rest.ensembl.org/"

class EnsemblRestQueries(object):
    '''Perform lookups using Ensembl's REST API'''

    def __init__(self, use_grch37_server=False, custom_server=None,
                 timeout=1.0, max_retries=2, reqs_per_sec=5,
                 log_level=logging.INFO):
        self._set_logger(logging_level=log_level)
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        if custom_server:
            self.server = custom_server
        elif use_grch37_server:
            self.server = grch37_server
        else:
            self.server = server
        self.timeout = timeout
        self.max_retries = max_retries

    def get_endpoint(self, endpoint, attempt=0):
        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                self.logger.debug("sleep {}s".format(1-delta))
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        self.logger.debug("Retrieving {}".format(self.server+endpoint))
        r = requests.get(self.server+endpoint, timeout=self.timeout,
                         headers={ "Content-Type" : "application/json"})
        self.req_count += 1
        if not r.ok:
            if attempt < self.max_retries:
                attempt += 1
                self.logger.info("Retry {}/{}".format(attempt,self.max_retries)
                                 + " for {}".format(self.server+endpoint))
                return self.get_endpoint(endpoint, attempt=attempt)
            r.raise_for_status()
        return r.json()

    def get_xref(self, query, all_levels='0', external_db=''):
        return self.get_endpoint("/xrefs/id/"+query+"?all_levels="+all_levels +
                                 ";external_db=" + external_db)

    def get_via_xref(self, query, species, get_type):
        endp = "/xrefs/symbol/{}/{}?object_type={}".format(species, query,
                                                           get_type)
        return self.get_endpoint(endp)

    def lookup_id(self, query, expand='0', phenotypes='0'):
        return self.get_endpoint("/lookup/id/"+query+"?expand="+expand+
                                 ";phenotypes="+phenotypes)

    def get_parent(self, query, expand='0'):
        info = self.lookup_id(query)
        if info:
            return self.lookup_id(info['Parent'], expand)
        return None

    def gene_from_enst(self, query, expand='0'):
        return self.get_parent(query, expand)

    def gene_from_ensp(self, query, expand='0'):
        trans = self.get_parent(query)
        if trans:
            return self.gene_from_enst(trans['id'])
        return None

    def lookup_variant(self, var, species='human', phenotypes='0'):
        return self.get_endpoint("/variation/"+species+"/"+var+"?phenotypes="
                                 + phenotypes)

    def lookup_ortholog(self, gene, taxon='10090'):
        ''' Lookup gene ID of ortholog (mouse by default) of given gene.'''
        data = self.get_endpoint("/homology/id/"+gene+"?target_taxon="+taxon)
        if data:
            try:
                return data['data'][0]['homologies'][0]['target']['id']
            except IndexError:
                return None
        return None

    def get_traits(self, gene):
        data = self.lookup_id(gene, phenotypes='1')
        if data and 'phenotypes' in data:
            return (phe['trait'] for phe in data['phenotypes'])
        return []

    def _set_logger(self, logging_level=logging.INFO):
        self.logger = logging.getLogger("Ensembl REST Queries")
        self.logger.setLevel(logging_level)
        formatter = logging.Formatter(
                        '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setLevel(self.logger.level)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

