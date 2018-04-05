'''
S = ABD.__str__()
s = re.search('(?<=e[\-\+]\d\d) ', S)

S_list = list(S)
while True:
    s = re.search('(?<=e[\-\+]\d\d) ', S)
    if not s: break
    S_list[s.start()] = ','
    S = ''.join(S_list)

while True:
    s = re.search('(?<=])\n', S)
    if not s: break
    a,b = s.span()
    S_list[a:b] = [',', '\n']
    print(s)
    S = ''.join(S_list)
    
while True:
    s = re.search('e[\-\+]00', S)
    if not s: break
    #S_list[s.start()] = ','
    #S = ''.join(S_list)
    a,b = s.span()
    S_list[a:b] = '    '
    S = ''.join(S_list)

S = S.replace('0.000', '0    ')

print(S)
'''
