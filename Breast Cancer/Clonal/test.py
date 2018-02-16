

if __name__ == '__main__':
    data  = [['B','1','15','0.8456'],['M','2','15','0.8456'],['M','3','15','0.8456'],['B','4','15','0.8456']]
    f_b = open('test_b','w')
    f_m = open('test_m','w')
    for element in data:
            if element[0] == 'B':
                phrase = ','.join(element)
                phrase += '\n'
                f_b.write(phrase)
    f_b.close()
    f_m.close()