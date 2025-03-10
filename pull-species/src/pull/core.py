from Bio import Entrez, SeqIO

from dotenv import load_dotenv

import os
import json

load_dotenv()
email = os.getenv("ENTREZ_EMAIL")
api_key = os.getenv("ENTREZ_API_KEY")

Entrez.email = email
Entrez.api_key = api_key

def search_species_sequences(species_name, gene_name="COX3", retmax=10):
  """
  Searches for nucleotide sequences of a specific gene in a given species.

  @param species_name (str): The scientific name of the species
  @param gene_name (str): The gene of interest (default is "COX3")
  @param retmax (int): Maximum number of records to retrieve (default is 10)

  @returns (list): A list of GenBank accession numbers
  """
  search_term = f"{species_name}[Organism] AND {gene_name}[Gene]"
  handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=retmax)
  record = Entrez.read(handle)
  handle.close()
  return record["IdList"]

def fetch_genbank(accession):
  """
  Fetches the GenBank record for a given accession number

  @param accession (str): The accession number of the desired GenBank record

  @returns (SeqRecord): A Biopython SeqRecord object containing the GenBank data
  """
  seq_record = None
  try:
    with Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text") as handle:
      seq_record = SeqIO.read(handle, "gb")
  except:
    print(f"Failed to fetch {accession}...")
  return seq_record

def cox3_translation_from_record(seq_record):
  """
  Grabs the translation of the COX3 gene from a GenBank SeqRecord, or None if
  not found...

  @param seq_record (SeqRecord): A Biopython SeqRecord object containing the 
    GenBank data

  @returns (str|None): The amino acid sequence of the COX3 gene translation, or
    None if not found
  """
  for feature in seq_record.features:
    if feature.type == "CDS" and "gene" in feature.qualifiers:
      if "COX3" in feature.qualifiers["gene"]:
        return feature.qualifiers.get("translation", ["<No translation available>"])[0]
  return None

def extract_cox3(species_name):
  """
  Grabs the translation of the COX3 gene given some species name...

  @param species_name (str): A latin species name to grab the COX3 gene sequence of

  @returns (str|None): The COX3 sequence for the provided species, or None if not
    found
  """

  gene_name = "COX3"
  accession_ids = search_species_sequences(species_name, gene_name)

  if not accession_ids or len(accession_ids) < 1:
    print(f"No {gene_name} gene sequences found for {species_name}.")
  else:
    accession_id = accession_ids[0]
    record = fetch_genbank(accession_id)
    if record is None:
      return None
    cox3_translation = cox3_translation_from_record(record)
    if cox3_translation is None:
      return None
    return cox3_translation

if __name__ == "__main__":
  species_list = [
    "Gorilla beringei",
    "Pan troglodytes",
    "Pongo pygmaeus",
    "Homo sapiens",
    "Lemur catta",
    "Panthera leo",
    "Panthera tigris",
    "Ursus arctos",
    "Canis lupus",
    "Felis catus",
    "Balaenoptera musculus",
    "Delphinus delphis",
    "Physeter macrocephalus",
    "Orcinus orca",
    "Tursiops truncatus",
    "Pteropus vampyrus",
    "Desmodus rotundus",
    "Myotis lucifugus",
    "Eonycteris spelaea",
    "Rhinolophus ferrumequinum",
    "Castor canadensis",
    "Hydrochoerus hydrochaeris",
    "Sciurus vulgaris",
    "Rattus norvegicus",
    "Cavia porcellus",
    "Equus ferus caballus",
    "Rhinoceros unicornis",
    "Tapirus terrestris",
    "Diceros bicornis",
    "Equus zebra",
    "Loxodonta africana",
    "Elephas maximus",
    "Trichechus manatus",
    "Dugong dugon",
    "Procavia capensis",
    "Orycteropus afer",
    "Dasypus novemcinctus",
    "Bradypus tridactylus",
    "Choloepus didactylus",
    "Myrmecophaga tridactyla",
    "Ornithorhynchus anatinus",
    "Tachyglossus aculeatus",
    "Macropus rufus",
    "Phascolarctos cinereus",
    "Vombatus ursinus",
    "Dendrolagus goodfellowi",
    "Sarcophilus harrisii",
    "Didelphis virginiana",
    "Monodelphis domestica",
    "Phascolosorex dorsalis",
    "Thylacinus cynocephalus",
    "Acinonyx jubatus",
    "Lynx lynx",
    "Puma concolor",
    "Leopardus pardalis",
    "Panthera onca",
    "Panthera uncia",
    "Mustela putorius furo",
    "Mephitis mephitis",
    "Procyon lotor",
    "Ailurus fulgens",
    "Enhydra lutris",
    "Odobenus rosmarus",
    "Mirounga leonina",
    "Halichoerus grypus",
    "Phoca vitulina",
    "Erinaceus europaeus",
    "Atelerix albiventris",
    "Sorex araneus",
    "Talpa europaea",
    "Condylura cristata",
    "Elephantulus rufescens",
    "Macroscelides proboscideus",
    "Solenodon paradoxus",
    "Tenrec ecaudatus",
    "Echinops telfairi",
    "Oryctolagus cuniculus",
    "Lepus europaeus",
    "Sylvilagus floridanus",
    "Ochotona princeps",
    "Camelus dromedarius",
    "Camelus bactrianus",
    "Vicugna vicugna",
    "Lama glama",
    "Bos taurus",
    "Bison bison",
    "Ovis aries",
    "Capra hircus",
    "Antilope cervicapra",
    "Gazella gazella",
    "Oryx dammah",
    "Alces alces",
    "Cervus elaphus",
    "Dama dama",
    "Giraffa camelopardalis",
    "Okapia johnstoni",
    "Hippopotamus amphibius",
    "Sus scrofa",
    "Phacochoerus africanus",
    "Dicotyles tajacu",
    "Bubalus bubalis",
    "Syncerus caffer",
    "Tragelaphus strepsiceros",
    "Taurotragus oryx",
    "Connochaetes taurinus",
    "Pantholops hodgsonii",
    "Rangifer tarandus",
    "Moschus moschiferus",
    "Capreolus capreolus",
    "Hydropotes inermis",
    "Neotragus pygmaeus",
    "Hippocamelus bisulcus",
    "Pudu puda",
    "Mazama americana",
    "Odocoileus virginianus",
    "Blastocerus dichotomus",
    "Ozotoceros bezoarticus",
    "Antilocapra americana",
    "Saiga tatarica",
    "Vicugna pacos",
    "Chinchilla lanigera",
    "Lagidium viscacia",
    "Erethizon dorsatum",
    "Coendou prehensilis",
    "Hydrochoerus isthmius",
    "Dasyprocta punctata",
    "Myocastor coypus",
    "Octodon degus",
    "Thryonomys swinderianus",
    "Petromus typicus",
    "Georychus capensis",
    "Heterocephalus glaber",
    "Spalax ehrenbergi",
    "Rattus rattus",
    "Mus musculus",
    "Peromyscus maniculatus",
    "Apodemus sylvaticus",
    "Clethrionomys glareolus",
    "Microtus arvalis",
    "Ondatra zibethicus",
    "Arvicola amphibius",
    "Lemmus lemmus",
    "Dicrostonyx torquatus",
    "Neofiber alleni",
    "Castor fiber",
    "Castor canadensis",
    "Hydromys chrysogaster",
    "Platypus australis",
    "Zaglossus bruijni",
    "Echidna hystrix",
    "Didelphis marsupialis",
    "Philander opossum",
    "Caluromys philander",
    "Monodelphis brevicaudata",
    "Thylamys elegans",
    "Dasyurus viverrinus",
    "Sminthopsis crassicaudata",
    "Antechinus stuartii",
    "Perameles nasuta",
    "Macrotis lagotis",
    "Notoryctes typhlops",
    "Vombatus hirsutus",
    "Lasiorhinus latifrons",
    "Phascolarctos cinereus",
    "Pseudocheirus peregrinus",
    "Petaurus breviceps",
    "Bettongia penicillata",
    "Aepyprymnus rufescens",
    "Dendrolagus matschiei",
    "Macropus giganteus",
    "Wallabia bicolor",
  ]

  # sequentially process each species to extract COX3 translation
  translations = {}

  for species in species_list:
    print(f"Fetching COX3 translation for {species}...")
    translation = extract_cox3(species)
    translations[species] = translation

  print("\n=== Summary of COX3 Translations ===")
  for species, translation in translations.items():
    if translation:
      print(f"[✓] {species}: Translation found, {len(translation)} amino acids")
    else:
      print(f"[✗] {species}: Translation not found")

  filtered_translations = {species: seq for species, seq in translations.items() if seq is not None}
  json_output = json.dumps(filtered_translations, indent=4)
  with open("cox3_translations.json", "w") as json_file:
    json_file.write(json_output)