import pandas as pd

df = pd.read_csv('ielex-subset.tsv', sep='\t', dtype=str)
df['concept'] = df['cc_alias'].apply(lambda s: s.split('-')[0])

out = ''


def clean_lexeme(lexeme):
	return lexeme.replace('  ', ' ').replace('  ', ' ').replace(' ', '/').lower()


def annotate_language(language, concept):
	x = df.loc[(df.language==language) & (df.concept==concept)]
	lexeme = x.lexeme.to_numpy()[0]
	lexeme = clean_lexeme(lexeme)
	return f'{lexeme}___{language}'


with open('wordtrees/ielex_contactrees.all.trees', 'r') as wordtree_file:
	line = wordtree_file.readline()
	while not line.strip() == 'Translate':
		out += line
		line = wordtree_file.readline()

	out += line
	line = wordtree_file.readline()
	
	while not line.strip() == ';':
		lang = line.split()[-1].strip(',')
		print((annotate_language(lang, 'all')))
		out += line.replace(lang, annotate_language(lang, 'all'))
		line = wordtree_file.readline()
		
	out += line
	for line in wordtree_file:
		out += line


with open('wordtrees_annotated/ielex_contactrees.all.trees', 'w') as wordtree_file:
	wordtree_file.write(out)
