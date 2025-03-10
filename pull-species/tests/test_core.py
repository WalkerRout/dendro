
from Bio.SeqFeature import SeqFeature, FeatureLocation

from unittest.mock import patch, MagicMock

import pytest

from pull.core import (
  search_species_sequences,
  fetch_genbank,
  cox3_translation_from_record,
  extract_cox3,
)

# mock
class FakeSeqRecord:
  def __init__(self, features):
    self.features = features

# helper for constructing 
def make_feature(gene_name="COX3", translation="MKT..."):
  return SeqFeature(
    type="CDS",
    location=FeatureLocation(0, 10),
    qualifiers={
      "gene": [gene_name],
      "translation": [translation],
    },
  )

def test_cox3_translation_from_record_found():
  record = FakeSeqRecord(features=[make_feature()])
  assert cox3_translation_from_record(record) == "MKT..."

def test_cox3_translation_from_record_not_found():
  record = FakeSeqRecord(features=[])
  assert cox3_translation_from_record(record) is None

@patch("pull.core.Entrez.esearch")
@patch("pull.core.Entrez.read")
def test_search_species_sequences(mock_read, mock_esearch):
  mock_read.return_value = {"IdList": ["ABC123", "XYZ789"]}
  mock_esearch.return_value = MagicMock()
  result = search_species_sequences("Homo sapiens", "COX3", 2)
  assert result == ["ABC123", "XYZ789"]
  mock_esearch.assert_called_once()

@patch("pull.core.SeqIO.read")
@patch("pull.core.Entrez.efetch")
def test_fetch_genbank_success(mock_efetch, mock_read):
  mock_handle = MagicMock()
  mock_efetch.return_value.__enter__.return_value = mock_handle
  mock_read.return_value = "FakeSeqRecord"
  result = fetch_genbank("ABC123")
  assert result == "FakeSeqRecord"

@patch("pull.core.SeqIO.read", side_effect=Exception("Fail"))
@patch("pull.core.Entrez.efetch")
def test_fetch_genbank_failure(mock_efetch, mock_read):
  mock_handle = MagicMock()
  mock_efetch.return_value.__enter__.return_value = mock_handle
  result = fetch_genbank("BADID")
  assert result is None

@patch("pull.core.fetch_genbank")
@patch("pull.core.search_species_sequences")
def test_extract_cox3_success(mock_search, mock_fetch):
  record = FakeSeqRecord(features=[make_feature("COX3", "MKT...")])
  mock_search.return_value = ["ABC123"]
  mock_fetch.return_value = record
  result = extract_cox3("Homo sapiens")
  assert result == "MKT..."

@patch("pull.core.search_species_sequences", return_value=[])
def test_extract_cox3_no_ids(mock_search):
  result = extract_cox3("Unknownus speciesus")
  assert result is None

@patch("pull.core.search_species_sequences", return_value=["BADID"])
@patch("pull.core.fetch_genbank", return_value=None)
def test_extract_cox3_no_record(mock_fetch, mock_search):
  result = extract_cox3("Homo sapiens")
  assert result is None

@patch("pull.core.search_species_sequences", return_value=["ABC123"])
@patch("pull.core.fetch_genbank")
def test_extract_cox3_no_translation(mock_fetch, mock_search):
  record = FakeSeqRecord(features=[make_feature("COX3", None)])
  mock_fetch.return_value = record
  record.features[0].qualifiers.pop("translation", None)
  result = extract_cox3("Homo sapiens")
  assert result == "<No translation available>"