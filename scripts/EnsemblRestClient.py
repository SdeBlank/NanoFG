import requests
import json
import time
import sys

class EnsemblRestClient(object):
    def __init__(self, reqs_per_sec=10):
        #self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        self.repeat = 0

    def perform_rest_action(self, server, endpoint, hdrs=None, parameters=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if parameters is None:
            parameters = {}
        data = None
        x=0
        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            sys.stderr.write("### RATE LIMIT - SELF")
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = requests.get(server + endpoint, headers=hdrs, params=parameters)
            request.raise_for_status()
            response = request.text
            if response:
                data = json.loads(response)
            self.req_count += 1

        except requests.exceptions.HTTPError as error:
            # check if we are being rate limited by the server
            if int(error.response.status_code) == 429:
                sys.stderr.write("### RATE LIMIT - SERVER")
                if 'Retry-After' in error.response.headers:
                    retry = error.response.headers['Retry-After']
                    time.sleep(float(retry))
                    data=self.perform_rest_action(server, endpoint, hdrs, parameters)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.response.status_code} Reason: {1.response.reason}\n'.format(server+endpoint, error))

        except requests.exceptions.ConnectionError as error:
            sys.stderr.write("### CONNECTION-ERROR")
            time.sleep(1)
            data=self.perform_rest_action(server, endpoint, hdrs, parameters)

        if data is None:
            sys.stderr.write("### DATA IS NONE - UNKNOWN ERROR")
            self.repeat += 1
            if self.repeat <= 5:
                time.sleep(4)
                data=self.perform_rest_action(server, endpoint, hdrs, parameters)
            else:
                sys.stderr.write("Too many tries to connect to the Ensembl database")
        return data
