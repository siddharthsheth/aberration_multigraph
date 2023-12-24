from aberration_multigraph.sv import BreakLocation, SVVertex, StructuralVariation, CNSegment
from collections import defaultdict
import os

class Patient:
    """
    A class to store a patient's data from the PCAWG project.
    This object stores both, copy number aberration data and structural
     variation data, for a patient.

    Attributes
    ----------
    id : str
        A patient's unique identifier string.
        Normally, a 128 bit string represented using 32 hexadecimal digits.
    dataset : str
        Determines the dataset containing the patient.
        Either 'tcga' or 'icgc'.
    logs : list
        A list of pairs of strings to store logging information and their
         priorities.
    bl_to_sv : defaultdict(dict)
        A dictionary to map BreakLocation to SVVertex.
        This prevents constructing multiple SVVertex objects for the same
         breakpoint location.
    sequence : defaultdict(list)
        A dictionary mapping chromosomes to lists of CNSegment objects that 
         store the sequencing information for this patient from the
         corresponding CNA file.
    svs : list
        A list of StructuralVariation objects to store all the SVs obeserved in
         this patient from the corresponding consensus SV file.
    bp_slack : int
        An additive error permitted when matching breakpoint locations in
         SVVertex and CNSegment.

    Methods
    -------
    _get_sequence(cn_path, cn_extension)
        Loads sequencing information for the patient.
    _get_svs(sv_path, sv_extension)
        Generates StructuralVariation objects for SVs observed in the patient.
    _chrom_to_int(chrom)
        Converts chromosome number from `str` to `int`.
    _create_sv_vertex(bl, pos, sv)
        Manage relationship between breakpoint locations and structural variations.
    check_valid()
        Checks if the data for this patient is clean.
    _check_sequence_complete()
        Checks if the sequencing data for this patient is complete.
    _bp_loc_cn(bl)
        Checks the copy number for the segment containing a break location.
    _bp_loc_rejoins(bl)
        Computes the number of rejoins a breakpoint location is involved in.
    _check_rejoins_cn()
        Checks if for every break location the number of rejoins exceeds the 
         number of copies.
    write_logs(log_file)
        Write logging info to file and clear logs.
    """

    def __init__(self, id, dataset, bp_slack=2):
        """
        Parameters
        ----------
        id : str
            A patient's unique identifier string.
            Normally, a 128 bit string represented using 32 hexadecimal digits.
        dataset : str
            Determines the dataset containing the patient.
            Either 'tcga' or 'icgc'.
        bp_slack : int, optional
            The acceptable absolute error when matching breakpoint locations
             from structural variations with copy number aberrations. 
        """
        sv_extension = '.pcawg_consensus_1.6.161116.somatic.sv.bedpe'
        cn_extension = '.consensus.20170119.somatic.cna.txt'
        sv_path = f'data/{dataset}_sv/open/'
        cn_path = f'data/{dataset}_cna/'
        
        self.id = id
        self.dataset = dataset
        self.logs = []
        self.bl_to_sv = defaultdict(dict)
        self.sequence = defaultdict(list)
        self.svs = []
        self.bp_slack = bp_slack

        self._get_svs(sv_path, sv_extension)
        self._get_sequence(cn_path, cn_extension)

    def _get_sequence(self, cn_path, cn_extension):
        """Loads sequencing information for the patient.
        
        This internal method uses the input parameters to open the CNA file for
         this patient and loads these sequences as `CNSegment` objects.
        Finally, all of these are stored as a list in the `sequence` attribute.
        This is the only time the method is called.

        This method raises a ValueError if the file cannot be opened or is in an
         incompatible format.

        Parameters
        ----------
        cn_path : str
            The relative path of the directory storing CNA file.
        cn_extension : str
            The extension of the CNA file (the part of the file name following
             the patient's ID.)

        Raises
        ------
        ValueError
            If the header does not match or CNA file cannot be opened.
        """
        file = self.id+cn_extension
        # sequence = defaultdict(list)
        with open(cn_path+file, 'r') as cnfile:                 #TODO: check if file can be opened
            headers = cnfile.readline()                         # read headers
            # if headers.split() != []                          #TODO: check headers
            for entry in cnfile:
                entry = entry.split()
                chrom = self._chrom_to_int(entry[0])
                for i, x in enumerate(entry[1:]):               # manage scientific notation
                    if 'e+' in x:
                        components = x.split('e+')
                        if len(components) == 2:
                            b, e = components
                            entry[i+1] = float(b)*10**int(e)
                segment = [x if x == 'NA' else int(x)  for x in entry[1:]]
                self.sequence[chrom].append(CNSegment(*segment))
        # return sequence
            
    def _get_svs(self, sv_path, sv_extension):
        """Generates StructuralVariation objects for SVs observed in the patient.
        
        This internal method uses the input parameters to open the consensus SV 
         file for this patient and loads these SVs as `StructuralVariation`
         objects.
        Finally, all of these are stored as a list in the `svs` attribute.
        This is the only time the method is called.

        Parameters
        ----------
        sv_path : str
            The relative path of the directory storing the consensus SV file.
        sv_extension : str
            The extension of the consensus SV file (the part of the file name
             following the patient's ID.)

        Raises
        ------
        ValueError
            If the header does not match or consensus SV file cannot be opened.
        """
        file = self.id+sv_extension
        # svs = []
        with open(sv_path+file, 'r') as svfile:                 #TODO: check if file can be opened
            headers = svfile.readline()                         # read headers
            # if headers.split() != []                          #TODO: check headers
            for entry in svfile:
                entry = entry.split()
                entry[0] = self._chrom_to_int(entry[0])
                entry[3] = self._chrom_to_int(entry[3])
                sv = StructuralVariation(*entry)
                self.svs.append(sv)
                bls = (BreakLocation(entry[0], int(entry[1])),
                        BreakLocation(entry[0], int(entry[2])),
                        BreakLocation(entry[3], int(entry[4])),
                        BreakLocation(entry[3], int(entry[5]))
                    )
                for i, bl in enumerate(bls):
                    self._create_sv_vertex(bl, i+1, sv)
        # return svs
    
    def _chrom_to_int(self, chrom):
        """Converts chromosome number from `str` to `int`.

        Parameters
        ----------
        chrom : str
            A string representing either a number from 1 to 22, or 'X', or 'Y'.

        Returns
        -------
        int
            An `int` representation of `chrom`.
            Inputs of 'X' or 'Y' are represented as 23 or 24 respectively.
        """
        if chrom == 'X':
            return 23
        elif chrom == 'Y':
            return 24
        else:
            return int(chrom)
        

    def _create_sv_vertex(self, bl, pos, sv):
        """Manage relationship between breakpoint locations and structural
         variations.

        If an SVVertex already exists for BreakLocation bl, the map stores
         StructuralVariation sv as another SV this SVVertex is involved in.
        Otherwise, a new SVVertex is created for bl and sv is stored as a
         StructuralVariation involving this BreakLocation.

        Parameters
        ----------
        bl : BreakLocation
            The breakpoint location (basis pair, chromosome) involved in a SV.
        pos : int
            The position of bl in the SV.
            Either 1, 2, 3, or 4.
        sv : StructuralVariation
            Object corresponding to the SV involving this break location.
        """
        # first check if SVVertex has already been created
        if bl in self.bl_to_sv:
            self.bl_to_sv[bl]['svs'].append((sv, pos))
        else:
            # create new SVVertex for this BreakLocation
            self.bl_to_sv[bl]['vertex'] = SVVertex(*bl)
            # update the StructuralVariation and the position in that SV this location belongs to
            self.bl_to_sv[bl]['svs'] = [(sv, pos)]


    def check_valid(self):
        """Checks if the data for this patient is clean.

        Any discrepancies observed while checking are stored in `logs`.
        """
        valid = True
        self.logs.append((5, f'Patient {self.id}\n'))

        if not self._check_sequence_complete():
            valid = False
        if not self._check_rejoins_cn():
            valid = False

        # for bl in self.bl_to_sv:
        #     print(f'{bl}: {self.bl_to_sv[bl]}')

        self.logs.append((5, '----------------------\n\n'))
        return valid

    def _check_sequence_complete(self):
        """Checks if the sequencing data for this patient is complete.

        It is acceptable if there are some missing base pairs at the beginning
         or the end of a chromosome.
        This method checks if there are any missing base pairs or duplicate base
         pairs in the sequenced data.
        
        Returns
        -------
        bool
            Returns True if sequence is complete, False otherwise.
        """
        old_logs = len(self.logs)
        for chrom in range(1,24):
            if chrom not in self.sequence:
                self.logs.append((5, f'Chromosome {chrom} not sequenced.\n'))
        for chrom in self.sequence:
            for i in range(1, len(self.sequence[chrom])):
                if self.sequence[chrom][i].start - self.sequence[chrom][i-1].end > 1:
                    self.logs.append((5, f"""Patient {self.id}: 
                                        Missing bps in chromosome {chrom}:
                                        {self.sequence[chrom][i-1].end} to
                                        {self.sequence[chrom][i].start}.\n"""))
                elif self.sequence[chrom][i].start - self.sequence[chrom][i-1].end < 1:
                    self.logs.append((5, f"""Patient {self.id}:
                                        Duplicate bps sequenced in chromosome {chrom}:
                                        {self.sequence[chrom][i].start} to
                                        {self.sequence[chrom][i-1].end}.\n"""))
        if len(self.logs) == old_logs:
            self.logs.append((5, f'Patient {self.id}: Sequencing data is complete.\n'))
            return True
        else:
            return False

    def _bp_loc_cn(self, bl):
        """Checks the copy number for the segment containing a break location.

        Given BreakLocation bl, this method returns the copy number of the
         CNSegment object containing that location

        Parameters
        ----------
        bl : BreakLocation
            The query breakpoint location

        Returns
        -------
        tuple
            The first element of the tuple is an integer indicating the total 
             copy number of the CNSegment containing the queried breakpoint.
            If the copy number was unavailable the first element is 0.
            The second element of the tuple is an error message.
            If the query was resolved satisfactorily, the error is 0.
            If the total copy number was marked 'NA', then the error is 'NA'.
            If the queried location was greater than the sequenced range, error
             is '>'.
            If the queried location was less than the sequenced range, error
             is '<'.
            Otherwise, the error is 'Not sequenced'.

        Raises
        ------
        ValueError
        TODO
        """
        chrom = bl.chrom
        bp = bl.bp
        for segment in self.sequence[chrom]:
            if segment.start <= bp <= segment.end:
                return (segment.total_cn, 0) if segment.total_cn != 'NA' else (0, 'NA')
        
        # log = f'Patient {self.id}: Copy number not found for break location {bl}.'
        if len(self.sequence[chrom]) == 0:
            logs.append((5, f'Unsequenced chromosome: {chrom} for {bl}'))
            return (0, 'Not sequenced')
        if bp < self.sequence[chrom][0].start:
            return (0, '<')
            # log += ' Location less than sequenced range.'
        elif bp > self.sequence[chrom][-1].end:
            return (0, '>')
            # log += ' Location greater than sequenced range.'
        # self.logs.append(log+'\n')
        # return (0, 'Not sequenced')
    
    def _bp_loc_rejoins(self, bl):
        """Computes the number of rejoins a breakpoint location is involved in.

        Given a breakpoint location, this method iterates over the
         StructuralVariation objects involving this location and returns the
         number of rejoins it was involved in.

        Parameters
        ----------
        bl : BreakLocation
            Represents a queried breakpoint location

        Returns
        -------
        int
            The number of rejoins this location was involved in.
        """
        rejoins = 0
        for sv, pos in self.bl_to_sv[bl]['svs']:
            if pos == 1 and sv.strand1 == '-':
                rejoins +=1
            elif pos == 2 and sv.strand1 == '+':
                rejoins += 1
            elif pos == 3 and sv.strand2 == '-':
                rejoins += 1
            elif pos == 4 and sv.strand2 == '+':
                rejoins += 1
        return rejoins

    def _check_rejoins_cn(self):
        """Checks if for every break location the number of rejoins exceeds the 
         number of copies.

        Returns
        -------
        bool
            True if for every break location the number of rejoins it is involved
             in does not exceed the number of copies of that location.
            False otherwise.
        """
        old_logs = len(self.logs)
        bp_locs = sorted([bl for bl in self.bl_to_sv])
        for bp_loc in bp_locs:
            rejoins = self._bp_loc_rejoins(bp_loc)
            copies, error = self._bp_loc_cn(bp_loc)
            if rejoins > copies:
                bp_less_slack = BreakLocation(bp_loc.chrom, bp_loc.bp-self.bp_slack)
                bp_less_slack_cn = self._bp_loc_cn(bp_less_slack)[0]
                bp_add_slack = BreakLocation(bp_loc.chrom, bp_loc.bp+self.bp_slack)
                bp_add_slack_cn = self._bp_loc_cn(bp_add_slack)[0]
                if error == 0:
                    if bp_less_slack_cn <= copies:
                        self.logs.append((1, f"""Patient {self.id}: {bp_loc} is rejoined
                                          in {rejoins} SVs while CN is {copies}.
                                          But a location -{self.bp_slack} is fine.\n"""))
                    elif  bp_add_slack_cn <= copies:
                        self.logs.append((1, f"""Patient {self.id}: {bp_loc} is rejoined
                                          in {rejoins} SVs while CN is {copies}.
                                          But a location +{self.bp_slack} is fine.\n"""))
                    else:
                        self.logs.append((5, f"""Patient {self.id}: {bp_loc} is rejoined
                                          in {rejoins} SVs while CN is {copies}.\n"""))
                elif error == '<':
                    if bp_add_slack_cn <= copies:
                        self.logs.append((1, f"""Patient {self.id}:
                                          {bp_loc} CN not found for break location.
                                          Location less than sequenced range.
                                          But a location +{self.bp_slack} is fine.\n"""))
                    else:
                        self.logs.append((5, f"""Patient {self.id}:
                                          {bp_loc} CN not found for break location.
                                          Location less than sequenced range.\n"""))
                elif error == '>':
                    if bp_less_slack_cn  <= copies:
                        self.logs.append((1, f"""Patient {self.id}:
                                          {bp_loc} CN not found for break location.
                                          Location greater than sequenced range.
                                          But a location -{self.bp_slack} is fine.\n"""))
                    else:
                        self.logs.append((5, f"""Patient {self.id}:
                                          {bp_loc} CN not found for break location.
                                          Location greater than sequenced range.\n"""))
                elif error == 'NA':
                    self.logs.append((5, f"""Patient {self.id}:
                                      {bp_loc} CN not found for break location.
                                      Marked NA.\n"""))
                else:
                    self.logs.append((5, f"""Patient {self.id}:
                                      {bp_loc} CN not found for break location.\n"""))
        if len(self.logs) == old_logs:
            self.logs.append((5, f"""Patient {self.id}: Success!
                              For every breakpoint #(rejoins) <= #(copies).\n"""))
            return True
        else:
            return False

    def write_logs(self, log_file, priority_level=5):
        """Write logging info to file and clear logs.

        Parameters
        ----------
        log_file : file
            The file to which the log data is appended.
        """
        with open(log_file, 'a') as log:
            log.writelines(log for priority, log in patient.logs
                            if priority >= priority_level)
        patient.logs = []


chromothripsis = '72f0a49a-aec8-47e5-846a-956c4da1507c.pcawg_consensus_1.6.161116.somatic.sv.bedpe'
simple = 'e1217ebe-1826-41a9-b6c4-702100a66f5e.pcawg_consensus_1.6.161116.somatic.sv.bedpe'
medium = '0ae2193f-0d68-485a-b8c2-7568cbcce33e.pcawg_consensus_1.6.161116.somatic.sv.bedpe'
files = [simple, medium]
# file += '.bedpe'
# path = '/Users/siddharthsheth/Dropbox/work/research-projects/TQFT Cancer Progression/tcga/open/'
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path+'/..')
# datasets = ['tcga', 'icgc']
datasets = ['tcga']
log_file = 'svlog.log'
if os.path.exists(log_file):
    os.remove(log_file)
logs = []
for dataset in datasets:
    sv_path = f'data/{dataset}_sv/open/'
    patient_ids = [file.split('.')[0] for file in os.listdir(sv_path)]
    patient_ids = [file.split('.')[0] for file in files]
    for i, id in enumerate(patient_ids):
        print(f'Working on patient number {i+1}: {id}.')
        patient = Patient(id, dataset)
        valid = patient.check_valid()
        # print(f'Working on patient number {i+1}: {id}. {valid}')
        patient.write_logs(log_file)
