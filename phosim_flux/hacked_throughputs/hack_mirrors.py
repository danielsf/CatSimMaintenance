for component in ('m1', 'm2', 'm3'):
    with open('%s.dat' % component, 'r') as in_file:
        in_lines = in_file.readlines()
    with open('%s.dat' % component, 'w') as out_file:
        for line in in_lines:
            if line[0] == '#':
                out_file.write(line)
                continue
            params = line.strip().split()
            out_file.write('%s 1.0\n' % params[0])
