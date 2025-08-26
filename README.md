Construção do modelo de baixa

1 - Ler a taxa de amostragem e o shape da sísmica
2 - Carregar os Horizontes, que devem estar no mesmo shape da sísmica
3 - Perfis em tempo
	Considere uma matriz com os poços empilhados um por linha contendo o mesmo tamanho e amostrado na mesma dimensão do modelo (nwells, nsamples)
	Normalmente esses dados ficam armazenados em arquivos ASCII ou .las. Será necessário ler e reamostrar para a mesma taxa da sísmica. Não pode conter NULL values
	Considere também uma matriz com duas colunas, cada uma com o index da localização de cada poço em il e xl, respectivamente. Não são os valores absolutos. Para construir essa matriz, será necessário extrair as localizações reais de cada poço e encontrar suas IL e XL correspondentes no dado sísmico. Ex: se IL for 5040, então é np.argmin(abs(IL - 5040))

4 - Inserir horizontes e construir o RGT (Relative Geological Time).
	Define-se o limite superior do modelo e o limite inferior
	O processo de criação do RGT é feito usando funções prontas e feito o seguinte
```python
mymodel = BMTools()
mymodel.build_rgt(surfaces=surfs, time=time, silence=False)
mymodel.rgt.shape
```

5 - Modelagem
    Definir um variograma
    Executar a interpolação (krigagem ordinária)
    Filtro passa baixa
    
Então os inputs são os horizontes (a superfície), perfis de poços (em duas matrizes que serão construídas com os poços em si, descrito no passo 3).
	
