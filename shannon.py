def get_frags(self, eta, iteration):
        if iteration == 0:
            for prim in self.prims:
                for frag in self.frags:
                    if prim in frag:
                        continue
                self.frags.append(list(prim))
        for i in range(0, len(self.frags)):
            self.frags[i] = sorted(set(self.frags[i]))
        oldfrags = copy.deepcopy(self.frags)
        for i in oldfrags:
            for atom in i:
                for j in self.prims:
                    for atom2 in j:
                        if self.A[atom][atom2] != 0:
                            if atom2 in self.frags:
                                continue
                            for atom3 in j:
                                for x in self.frags:
                                    x.append(atom3)
        for i in range(0, len(self.frags)):
            self.frags[i] = set(self.frags[i])

        iteration = iteration + 1
        if eta > iteration:
            print('d')
            self.get_frags(eta, iteration)


