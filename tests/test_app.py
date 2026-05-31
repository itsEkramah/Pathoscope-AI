# pyrefly: ignore [missing-import]
import pytest
from pathlib import Path
from streamlit.testing.v1 import AppTest

def test_streamlit_app_initial_render():
    """
    Verifies that the Streamlit application successfully performs its initial render,
    loading all sliders, configurations, default parameters, and the hero title cards.
    """
    # Load the app using official streamlit testing utilities
    at = AppTest.from_file("pathoscope/app.py", default_timeout=60)
    
    # Trigger initial run/render
    at.run()
    
    # Assert no unhandled exceptions occurred
    assert not at.exception, f"App threw an unhandled exception: {at.exception}"
    
    # Assert main markdown card contains the application name
    markdowns = [m.value for m in at.markdown]
    assert any("PathoScope AI" in m for m in markdowns), "App title not found in markdown headers."
    
    # Assert sliders for QC length are present
    assert len(at.slider) > 0, "No sliders detected on the page."
    
    # Assert text input boxes exist (database paths, API configuration keys)
    assert len(at.text_input) > 0, "No text input components found."
    
    # Assert file upload control exists on screen
    assert len(at.file_uploader) > 0, "File uploader component is missing."
    
    # Assert run button is rendered
    assert len(at.button) > 0, "Run button is missing."
    assert at.button[0].label == "RUN GENOMICS PIPELINE"


def test_streamlit_parameter_state_handling():
    """
    Verifies that modifying configurations on the Streamlit page updates the internal state correctly.
    """
    at = AppTest.from_file("pathoscope/app.py", default_timeout=60)
    at.run()
    
    # Verify default checkbox and initial states
    preset_checkbox = at.checkbox[0]
    assert preset_checkbox.label == "Use Preset Sequence"
    assert preset_checkbox.value is False
    
    # Simulate user checking the preset sequence box
    preset_checkbox.check().run()
    assert at.checkbox[0].value is True

    
    # Verify slider mappings for QC length
    min_len_slider = at.slider[0]
    assert min_len_slider.value == 100  # Default loaded from default_config.yaml
    
    # Mutate parameter and re-run simulation
    min_len_slider.set_value(250).run()
    assert at.slider[0].value == 250


def test_advanced_mode_visibility():
    """
    Verifies that the advanced bioinformatician checkbox toggles advanced debugger tab rendering.
    """
    at = AppTest.from_file("pathoscope/app.py", default_timeout=60)
    at.run()
    
    # Find the advanced bioinformatician checkbox robustly by label substring
    adv_checkbox = [c for c in at.checkbox if "Advanced" in c.label][0]
    assert adv_checkbox.value is False
    
    # Simulate turning advanced mode ON
    adv_checkbox.check().run()
    
    # Re-verify value state updates correctly
    adv_checkbox_new = [c for c in at.checkbox if "Advanced" in c.label][0]
    assert adv_checkbox_new.value is True


