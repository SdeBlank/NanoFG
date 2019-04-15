import requests
import json
import time
import sys

class EnsemblRestClient(object):
    def __init__(self, reqs_per_sec=15):
        #self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        self.repeat = 0
        self.req_times=[]

    def perform_rest_action(self, server, endpoint, hdrs=None, parameters=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if parameters is None:
            parameters = {}
        data = None
        # check if we need to rate limit ourselves
        if len(self.req_times)==self.reqs_per_sec:
            print(abs(self.req_times[0]-self.req_times[self.reqs_per_sec-1]))
            if abs(self.req_times[0]-self.req_times[self.reqs_per_sec-1]) < 1:
                delta = time.time() - self.last_req
                if delta < 1:
                    time.sleep(1 - delta)
                self.last_req = time.time()
                self.req_count = 0

        try:
            if len(self.req_times)>=self.reqs_per_sec:
                del self.req_times[0]
            self.req_times.append(time.time())
            request = requests.get(server + endpoint, headers=hdrs, params=parameters)
            request.raise_for_status()
            response = request.text
            if response:
                data = json.loads(response)
            self.req_count += 1

        except requests.exceptions.ConnectionError as error:
            time.sleep(1)
            data=self.perform_rest_action(server, endpoint, hdrs, parameters)

        except requests.exceptions.HTTPError as error:
            # check if we are being rate limited by the server
            if int(error.response.status_code) == 429:
                if 'Retry-After' in error.response.headers:
                    retry = error.response.headers['Retry-After']
                    time.sleep(float(retry))
                    data=self.perform_rest_action(server, endpoint, hdrs, parameters)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.response.status_code} Reason: {1.response.reason}\n'.format(server+endpoint, error))

        if data is None:
            self.repeat += 1
            if self.repeat <= 5:
                time.sleep(4)
                data=self.perform_rest_action(server, endpoint, hdrs, parameters)
            else:
                sys.stderr.write("Too many tries to connect to the Ensembl database\n")
        return data
