

CELTIC = [
    'Old_Irish', 'Irish_A', 'Scots_Gaelic', 'Welsh_N', 'Welsh_C', 'Breton_List', 'Breton_Se',
    'Breton_ST', 'Cornish', 'Irish_B', 'Manx', 'Old_Breton', 'Old_Welsh', 'Gaulish', 'Old_Cornish'
]
GERMANIC = [
    'Gothic', 'Swedish_Up', 'Swedish_Vl', 'Swedish', 'Danish', 'Old_Norse', 'Icelandic_ST', 'Faroese', 'Dutch_List',
    'Flemish', 'Frisian', 'Old_High_German', 'German', 'Luxembourgish', 'Old_English', 'English', 'German_Munich',
    'Afrikaans', 'Norwegian', 'Pennsylvania_Dutch', 'Schwyzerdutsch', 'Elfdalian'
]
ROMANCE = [
    'Umbrian', 'Oscan',  'Latin', 'Provencal', 'French', 'Walloon', 'Spanish', 'Portuguese_ST', 'Catalan', 'Sardinian_N',
    'Sardinian_L', 'Sardinian_C', 'Rumanian_List', 'Vlach', 'Italian', 'Ladin', 'Friulian', 'Romansh', 'Dolomite_Ladino'
]

BALTO_SLAVIC = [
    'Old_Prussian', 'Lithuanian_ST', 'Latvian', 'Slovenian', 'Lower_Sorbian', 'Upper_Sorbian', 'Czech', 'Slovak',
    'Czech_E', 'Ukrainian', 'Byelorussian', 'Polish', 'Russian', 'Macedonian', 'Bulgarian', 'Serbocroatian',
    'Old_Church_Slavonic', 'BULGARIAN_P', 'BYELORUSSIAN_P', 'CZECH_P', 'Lithuanian_O', 'MACEDONIAN_P',
    'RUSSIAN_P', 'POLISH_P', 'Serbian', 'SERBOCROATIAN_P', 'SLOVAK_P', 'SLOVENIAN_P', 'UKRAINIAN_P'
]
ALBANIAN = [
    'ALBANIAN', 'Albanian_C', 'Albanian_G', 'Albanian_Standard', 'Albanian_T', 'Albanian_Top', 'Albanian_K'
]
ARMENIAN = [
    'Armenian_List', 'Armenian_Mod', 'Classical_Armenian'
]
GREEK = [
    'Greek_D', 'Greek_K', 'Greek_Md', 'Greek_Ml', 'Greek_Mod', 'Ancient_Greek'
]
IRANIAN = [
    'Baluchi', 'Avestan', 'Kurdish', 'Persian', 'Tadzik', 'Zazaki', 'Digor_Ossetic',
    'Iron_Ossetic', 'Sogdian', 'Wakhi', 'Pashto', 'Shughni', 'Sariqoli', 'Ossetic'
]
INDO_ARYAN = [
    'Kashmiri', 'Bengali', 'Bihari', 'Hindi', 'Khaskura', 'Lahnda', 'Marathi', 'Marwari',
    'Nepali', 'Oriya', 'Panjabi_ST', 'Sindhi', 'Urdu', 'Old_Persian', 'Vedic_Sanskrit',
    'Waziri', 'Assamese', 'Gujarati', 'Magahi', 'Singhalese'
]
OTHER = [
    'Tocharian_A', 'Tocharian_B', 'Hittite', 'Luvian', 'Lycian', 'Palaic', 'Kati'
]

EXCLUDE = [
    'Sranan', 'Brazilian', 'French_Creole_C', 'French_Creole_D', 'Gypsy_Gk', 'Proto-Indo-European'
]


CHANG_NARROW = []

CHANG_MEDIUM = [
    'Afrikaans', 'Albanian_K', 'Albanian_Top', 'Ancient_Greek', 'Armenian_List',
    'Armenian_Mod', 'Assamese', 'Avestan', 'Baluchi', 'Bengali', 'Bihari', 'Breton_ST',
    'Bulgarian', 'Byelorussian', 'Catalan', 'Classical_Armenian', 'Cornish', 'Czech',
    'Danish', 'Digor_Ossetic', 'Dutch_List', 'English', 'Faroese', 'Flemish', 'French',
    'Frisian', 'Friulian', 'German', 'Gothic', 'Greek_Mod', 'Gujarati', 'Gypsy_Gk',
    'Hindi', 'Hittite', 'Icelandic_ST', 'Irish_B', 'Italian', 'Kashmiri', 'Ladin',
    'Lahnda', 'Latin', 'Latvian', 'Lithuanian_ST', 'Luxembourgish', 'Macedonian',
    'Marathi', 'Nepali', 'Norwegian', 'Old_Church_Slavonic', 'Old_English',
    'Old_High_German', 'Old_Irish', 'Old_Norse', 'Oriya', 'Panjabi_ST', 'Pashto',
    'Persian', 'Polish', 'Portuguese_ST', 'Provencal', 'Romansh', 'Rumanian_List',
    'Russian', 'Sardinian_C', 'Sardinian_N', 'Schwyzerdutsch', 'Scots_Gaelic', 'Serbian',
    'Singhalese', 'Slovak', 'Slovenian', 'Spanish', 'Swedish', 'Tadzik', 'Tocharian_B',
    'Ukrainian', 'Upper_Sorbian', 'Urdu', 'Vedic_Sanskrit', 'Walloon', 'Waziri', 'Welsh_N']

CHANG_BROAD = []




# RENAMING DICTS (necessary to translate between different naming conventions)
RENAME_CHANG = {
    'Adapazar': 'Armenian_List',
    'Arvanitika': 'Albanian_K',
    'Belarusian': 'Byelorussian',
    'Breton': 'Breton_ST',
    'Cagliari': 'Sardinian_C',
    'Dutch': 'Dutch_List',
    'Eastern_Armenian': 'Armenian_Mod',
    'Icelandic': 'Icelandic_ST',
    'Irish': 'Irish_B',
    'Lithuanian': 'Lithuanian_ST',
    'Modern_Greek': 'Greek_Mod',
    'Nuorese': 'Sardinian_N',
    'Old_Church_Slavic': 'Old_Church_Slavonic',
    'Old_Irish_B': 'Old_Irish',
    'Old_West_Norse': 'Old_Norse',
    'Panjabi': 'Panjabi_ST',
    'Portuguese': 'Portuguese_ST',
    'Romani': 'Gypsy_Gk',
    'Romanian': 'Rumanian_List',
    'Swiss_German': 'Schwyzerdutsch',
    'Tajik': 'Tadzik',
    'Tosk': 'Albanian_Top',
    'Welsh': 'Welsh_N',
}

RENAME_BOUCKAERT = {
    'Breton_Se': 'Breton_SE',
    'Rumanian_List': 'Romanian_List',
    'Swedish_Vl': 'Swedish_VL',
    'Swedish': 'Swedish_List',
    'German': 'German_ST',
    'English': 'English_ST',
    'Lower_Sorbian': 'Lusatian_L',
    'Upper_Sorbian': 'Lusatian_U',
    'Old_Church_Slavonic': 'Old_Church_Slavic',
    'Greek_Ml': 'Greek_ML',
}
