from abc import ABC, abstractmethod


class BiologicalSequence(ABC):
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __getitem__(self, item):
        pass

    @abstractmethod
    def alphabet_check(self):
        pass


class Sequence(BiologicalSequence):
    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return str(self.sequence)

    def __getitem__(self, slc):
        return self.sequence[slc]

    def alphabet_check(self):
        alphabet = set(self.sequence)
        if not alphabet.issubset(self.allowed):
            raise ValueError('Wrong alphabet')


class NucleicAcidSequence(Sequence):
    def __init__(self, sequence):
        self.sequence = sequence
        self.allowed = {'A', 'T', 'G', 'C', 'U'}

    def complement(self):
        self.complemented = []
        for letter in self.sequence:
            self.complemented.append(self.complement_dict[letter])
        return type(self)(''.join(self.complemented))

    def gc_content(self):
        g_count = self.sequence.count('G')
        c_count = self.sequence.count('C')
        content = (g_count + c_count) / self.sequence.__len__()
        return content


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence):
        super().__init__(sequence)
        self.allowed = {'A', 'T', 'G', 'C'}
        self.complement_dict = {'A': 'T', 'T': 'A',
                                'G': 'C', 'C': 'G'}

    def transcribe(self):
        return self.sequence.replace('T', 'U')


class RNASequence(NucleicAcidSequence):
    def __init__(self, sequence):
        super.__init__(sequence)
        self.allowed = {'A', 'C', 'G', 'U'}
        self.complement_dict = {'A': 'U', 'U': 'A',
                                'G': 'C', 'C': 'G'}


class AminoAcidSequence(Sequence):
    def __init__(self, sequence):
        self.sequence = sequence
        self.allowed = {'A', 'R', 'N', 'D',
                        'C', 'Q', 'E', 'G',
                        'H', 'I', 'L', 'K',
                        'M', 'F', 'P', 'S',
                        'T', 'W', 'Y', 'V'}

        self.composition = {'A': 0, 'R': 0, 'N': 0, 'D': 0,
                            'C': 0, 'Q': 0, 'E': 0, 'G': 0,
                            'H': 0, 'I': 0, 'L': 0, 'K': 0,
                            'M': 0, 'F': 0, 'P': 0, 'S': 0,
                            'T': 0, 'W': 0, 'Y': 0, 'V': 0}

    def count_composition(self):
        for am_acid in self.sequence:
            self.composition[am_acid] += 1
        for am_acid in self.composition:
            self.composition[am_acid] = round(self.composition[am_acid] * 100 / len(self.sequence), 2)

        return self.composition
