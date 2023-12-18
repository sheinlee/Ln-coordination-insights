import os
element = 'La'

os.system("python ./"+ element +"_phenanthroline.py")
os.system("python ./"+ element +"_phenanthroline2.py")
os.system("python ./"+ element +"_phenanthroline3.py")
os.system("python ./"+ element +"_phenanthroline4.py")



f1 = open(element+'_phenanthroline_OK.txt','r')
f2 = open(element+'_phenanthroline2_OK.txt','r')
f3 = open(element+'_phenanthroline3_OK.txt','r')
f4 = open(element+'_phenanthroline4_OK.txt','r')

txt1 = f1.readlines()
txt2 = f2.readlines()
txt3 = f3.readlines()
txt4 = f4.readlines()

with open(element+'_phenanthroline_OK_1.txt','a') as f:
    f.close()
with open(element+'_phenanthroline_OK_2.txt','a') as f:
    f.close()
with open(element+'_phenanthroline_OK_3.txt','a') as f:
    f.close()
with open(element+'_phenanthroline_OK_4.txt','a') as f:
    f.close()

for w in txt1:
    with open(element+'_phenanthroline_OK_1.txt','a') as f:
        if w not in txt2:
            f.write(w)
        else:
            print("已去除重复--> "+w)
        f.close()
f1.close()  

for w in txt2:
    with open(element+'_phenanthroline_OK_2.txt','a') as f:
        if w not in txt3:
            f.write(w)
        else:
            print("已去除重复--> "+w)
        f.close()
f2.close()  

for w in txt3:
    with open(element+'_phenanthroline_OK_3.txt','a') as f:
        if w not in txt4:
            f.write(w)
        else:
            print("已去除重复--> "+w)
        f.close()
f3.close()  

f4.close()


     

#将.txt文件转化成.gcd文件
os.rename(element+'_phenanthroline_OK_1.txt',element+'_phenanthroline_OK_1.gcd')
os.rename(element+'_phenanthroline_OK_2.txt',element+'_phenanthroline_OK_2.gcd')
os.rename(element+'_phenanthroline_OK_3.txt',element+'_phenanthroline_OK_3.gcd')
os.rename(element+'_phenanthroline4_OK.txt',element+'_phenanthroline_OK_4.gcd')