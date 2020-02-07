'''
Copyright © 2018 Anton Tsukanov. Contacts: tsukanov@bionet.nsc.ru
License: http://www.gnu.org/licenses/gpl.txt

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

    ########################
    #FIND MODEL BY CHIPMUNK#
    ########################
    if not os.path.isfile(chipmunk + '/CHIPMUNK_MOTIF.txt'):
        #Create fastaWig for chipmunk
        print('chipmunk find motifs for {0}'.format(tag))
        run_chipmunk_fasta(path_to_java, path_to_chipmunk,
        fasta + '/' + tag + '_'+ str(training_sample_size) + '.fa',
        chipmunk + '/CHIPMUNK_MOTIF.txt',
        motif_length_start, motif_length_end,
        try_size, cpu_count, zoops)
    else:
        print('File {0} already exists'.format(chipmunk + '/chipmunk_MOTIF.txt'))

    ###########################################################################
    #Parse results of chipmunk into files .meme, .pwm and .fasta (multi fasta)#
    ###########################################################################

    args = ['python3', path_to_python_tools + 'parse_chipmunk_results.py',
           '-i', chipmunk + '/CHIPMUNK_MOTIF.txt',
           '-o', chipmunk,
           '-t', tag + '_' + 'CHIPMUNK_MOTIF']
    r = subprocess.call(args)


    ##############################################################################
    #Get oPWM from chipmunk results. OUTPUT: .meme, .pwm and .fasta (multi fasta)#
    ##############################################################################
    if not os.path.isfile(motifs + '/' + tag + '_' + 'OPTIMAL_MOTIF.meme'):
        args = ['python3', path_to_python_tools + 'make_oPWM.py',
                '-c', chipmunk + '/CHIPMUNK_MOTIF.txt',
                '-f', fasta + '/' + tag + '_'+ str(training_sample_size) + '.fa',
                '-n', '5000',
                '-P', cpu_count,
                '-o', motifs,
                '-t', tag + '_' + 'OPTIMAL_MOTIF']
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(motifs + '/PEAKS039334_OPTIMAL_MOTIF.meme'))


    ##################
    #GET MOTIF LENGTH#
    ##################

    with open(motifs + '/' + tag + '_' + 'OPTIMAL_MOTIF.fasta', 'r') as file:
        for i in file:
            if i.startswith('>'):
                continue
            else:
                motif_length = len(i.strip())
                break
    file.close()