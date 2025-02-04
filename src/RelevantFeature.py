"""Class to be used with NISTSearch class where the search results will be stored as a dictionary of RelevantFeature
    objects to later be used with the IntensityCollector class to format output and determine what data to collect"""

class RelevantFeature:

    def __init__(self, name: str, retention_time: float, top_ions: list, match_score: float):
        self.name = name
        self.retention_time = retention_time
        self.top_ions = top_ions
        self.match_score = match_score

    def __repr__(self):
        return f"RelevantFeature(name = {self.name}, retention_time = {self.retention_time}, top_ions = {self.top_ions}, match_score = {self.match_score})"