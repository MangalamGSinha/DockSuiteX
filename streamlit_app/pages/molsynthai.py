"""
DockSuiteX — MolSynthAI
Embeds the MolSynthAI lead-optimization tool directly inside DockSuiteX.
"""
import streamlit as st


# ── Style: hide outer scrollbar, make iframe fill viewport ──
st.markdown(
    """
    <style>
        /* Hide the outer Streamlit scrollbar */
        .stMainBlockContainer, .main .block-container,
        section.main, section.main > div {
            overflow: hidden !important;
        }
        /* Remove padding & max-width */
        .block-container {
            padding: 0 !important;
            max-width: 100% !important;
        }
        /* Borderless, full-width iframe */
        iframe {
            border: none !important;
            width: 100% !important;
        }
    </style>
    """,
    unsafe_allow_html=True,
)

# ── Embed MolSynthAI (use viewport-safe height) ──
st.iframe(
    src="https://molsynthai.streamlit.app/?embedded=true",
    height=800,
)

# Inject JS: resize iframe to fill viewport & ensure it scrolls internally
st.html(
    """
    <script>
    (function() {
        var frames = window.parent.document.querySelectorAll('iframe');
        for (var i = 0; i < frames.length; i++) {
            var src = frames[i].getAttribute('src') || '';
            if (src.indexOf('molsynthai') !== -1) {
                frames[i].setAttribute('scrolling', 'yes');
                frames[i].style.height = 'calc(100vh - 3.5rem)';
                frames[i].style.minHeight = '600px';
                break;
            }
        }
    })();
    </script>
    """,
)
